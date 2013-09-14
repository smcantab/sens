
from sens import _database_normalmodes
import sqlalchemy
import sys

import argparse

def add_normalmodes_table(connection):
    ''' migrating from version 0 to 1
    
        fields for log product of frequencies and point group order were added
        to Minimum and TransitionState
    '''
    print "adding normalmodes table"
    connection.execute("""CREATE TABLE tbl_normal_modes (
    freqs BLOB,
    vectors BLOB,
    _minimum_id INTEGER NOT NULL,
    PRIMARY KEY (_minimum_id),
    FOREIGN KEY(_minimum_id) REFERENCES tbl_minima (_id)
    );""")


def migrate(db):
    engine = sqlalchemy.create_engine("sqlite:///%s"%db)
    
    
    connection = engine.connect()
    trans = connection.begin()
    try:
        add_normalmodes_table(connection)
    except RuntimeError:
        trans.rollback()
        raise#raise RuntimeError("failed to migrate database")
    trans.commit()
    connection.close()

def main():    
    parser = argparse.ArgumentParser(description="add the NormalModes table to an existing database")
    parser.add_argument("dbfile", type=str, help="database file name")
    args = parser.parse_args()

    migrate(args.dbfile)

if __name__ == '__main__':
    main()