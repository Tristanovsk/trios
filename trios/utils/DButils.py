import psycopg2

username, pwd = 'tristan.harmel', '7.zcBjFa' #*****'

IDlake = 'SRC04'

DBserver = 'serveurbd.aix.irstea.priv'
DBname = 'bd_plando'

db = psycopg2.connect(
    "host= " + DBserver + " port='5434' dbname='" + DBname +
    "' user='" + username + "' password='" + pwd + "'")

SQL = "SELECT altitude_pla FROM plan_eau WHERE code_lac='" + IDlake + "';"

altitude = pgSQLquery2numpy(db, SQL)['altitude_pla'][0]
