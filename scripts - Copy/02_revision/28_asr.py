import json
sites=json.load(open('/home/claude/pet/analysis/catalytic_sites.json'))
print(json.dumps(sites,indent=1)[:1200])
