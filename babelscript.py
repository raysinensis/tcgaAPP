pybabel extract -F babel.cfg -o messages.pot .
pybabel init -i messages.pot -d translations -l zh_Hans_CN
#pybabel update -i messages.pot -d translations
pybabel compile -d translations
