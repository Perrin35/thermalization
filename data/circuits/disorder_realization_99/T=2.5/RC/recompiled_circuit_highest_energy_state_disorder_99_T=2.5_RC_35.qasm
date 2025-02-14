OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9361629) q[0];
sx q[0];
rz(-0.054447629) q[0];
sx q[0];
rz(0.86583889) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39393878) q[0];
sx q[0];
rz(-2.9705715) q[0];
sx q[0];
rz(-1.8500516) q[0];
rz(0.073926386) q[2];
sx q[2];
rz(-1.0763028) q[2];
sx q[2];
rz(-1.489515) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9283951) q[1];
sx q[1];
rz(-1.0093912) q[1];
sx q[1];
rz(2.2092845) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9179814) q[3];
sx q[3];
rz(-1.2652186) q[3];
sx q[3];
rz(-1.7479727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8286579) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(0.9642967) q[2];
rz(-3.0621081) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(-2.3704119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013414772) q[0];
sx q[0];
rz(-2.7981813) q[0];
sx q[0];
rz(-1.5403904) q[0];
rz(0.84114289) q[1];
sx q[1];
rz(-1.7336188) q[1];
sx q[1];
rz(-0.85743633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4682789) q[0];
sx q[0];
rz(-0.43860093) q[0];
sx q[0];
rz(-1.5836141) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.964999) q[2];
sx q[2];
rz(-0.73179663) q[2];
sx q[2];
rz(0.082187637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.407436) q[1];
sx q[1];
rz(-0.70824558) q[1];
sx q[1];
rz(0.70869653) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8624865) q[3];
sx q[3];
rz(-1.4785023) q[3];
sx q[3];
rz(-1.5444322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2555344) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(1.9319755) q[2];
rz(-1.7602734) q[3];
sx q[3];
rz(-1.2608903) q[3];
sx q[3];
rz(2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7968314) q[0];
sx q[0];
rz(-1.3153356) q[0];
sx q[0];
rz(-1.8094081) q[0];
rz(1.6650797) q[1];
sx q[1];
rz(-1.3308728) q[1];
sx q[1];
rz(0.61980334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22794321) q[0];
sx q[0];
rz(-1.648128) q[0];
sx q[0];
rz(-1.5864775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1479883) q[2];
sx q[2];
rz(-2.8393778) q[2];
sx q[2];
rz(-0.72590088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9732144) q[1];
sx q[1];
rz(-1.4992979) q[1];
sx q[1];
rz(2.8236103) q[1];
x q[2];
rz(1.2215106) q[3];
sx q[3];
rz(-2.4716645) q[3];
sx q[3];
rz(1.0494029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1406113) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(-0.78835431) q[2];
rz(-1.8654478) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(-2.2787794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514423) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(2.0748806) q[0];
rz(-1.7881296) q[1];
sx q[1];
rz(-2.2746494) q[1];
sx q[1];
rz(1.9852759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8114093) q[0];
sx q[0];
rz(-1.5519987) q[0];
sx q[0];
rz(1.4292595) q[0];
x q[1];
rz(-0.026211003) q[2];
sx q[2];
rz(-1.6115341) q[2];
sx q[2];
rz(-2.7762976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51976065) q[1];
sx q[1];
rz(-1.9433738) q[1];
sx q[1];
rz(-0.31078679) q[1];
x q[2];
rz(1.636999) q[3];
sx q[3];
rz(-0.90747031) q[3];
sx q[3];
rz(2.5732793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(0.058825292) q[2];
rz(-0.63932738) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-2.2842469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0462129) q[0];
sx q[0];
rz(-1.7261427) q[0];
sx q[0];
rz(3.131102) q[0];
rz(1.563021) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(-0.48935997) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45324126) q[0];
sx q[0];
rz(-2.4972389) q[0];
sx q[0];
rz(-2.8768566) q[0];
rz(-pi) q[1];
rz(2.0496558) q[2];
sx q[2];
rz(-0.74981028) q[2];
sx q[2];
rz(0.15287003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4976144) q[1];
sx q[1];
rz(-2.3757739) q[1];
sx q[1];
rz(0.52180565) q[1];
rz(0.030767083) q[3];
sx q[3];
rz(-1.5813785) q[3];
sx q[3];
rz(-2.2969674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2995149) q[2];
sx q[2];
rz(-2.044951) q[2];
sx q[2];
rz(-0.00046029885) q[2];
rz(-0.67356235) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(-2.4323997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29086581) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(-1.3036183) q[0];
rz(0.29779008) q[1];
sx q[1];
rz(-2.2460263) q[1];
sx q[1];
rz(-2.8477125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0336908) q[0];
sx q[0];
rz(-1.0009888) q[0];
sx q[0];
rz(0.66500647) q[0];
x q[1];
rz(1.05386) q[2];
sx q[2];
rz(-0.74956761) q[2];
sx q[2];
rz(0.38656879) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0386427) q[1];
sx q[1];
rz(-1.0979196) q[1];
sx q[1];
rz(2.6421618) q[1];
rz(-0.31556337) q[3];
sx q[3];
rz(-1.7949008) q[3];
sx q[3];
rz(-1.4423646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46875724) q[2];
sx q[2];
rz(-1.4361278) q[2];
sx q[2];
rz(-0.22077665) q[2];
rz(1.2729493) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(1.1481736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64977589) q[0];
sx q[0];
rz(-0.62115541) q[0];
sx q[0];
rz(2.5671) q[0];
rz(2.5993787) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(2.7508459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3706659) q[0];
sx q[0];
rz(-0.61425137) q[0];
sx q[0];
rz(0.25247981) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0263881) q[2];
sx q[2];
rz(-1.6154535) q[2];
sx q[2];
rz(-2.3400729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73394164) q[1];
sx q[1];
rz(-1.9875257) q[1];
sx q[1];
rz(1.5647566) q[1];
rz(-pi) q[2];
rz(1.5338906) q[3];
sx q[3];
rz(-1.9362984) q[3];
sx q[3];
rz(2.9707462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6812402) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(-2.528842) q[2];
rz(-0.98226205) q[3];
sx q[3];
rz(-1.3145612) q[3];
sx q[3];
rz(-2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9923582) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(3.0862578) q[0];
rz(-1.5570359) q[1];
sx q[1];
rz(-2.0535856) q[1];
sx q[1];
rz(2.7016644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81145168) q[0];
sx q[0];
rz(-2.6127671) q[0];
sx q[0];
rz(-2.6374966) q[0];
rz(-pi) q[1];
rz(1.7753168) q[2];
sx q[2];
rz(-2.0530982) q[2];
sx q[2];
rz(-2.4647146) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7288558) q[1];
sx q[1];
rz(-1.8860894) q[1];
sx q[1];
rz(-0.58087491) q[1];
rz(-pi) q[2];
rz(-1.004435) q[3];
sx q[3];
rz(-1.8105339) q[3];
sx q[3];
rz(-0.15739239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3736734) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(-2.8295753) q[2];
rz(-0.9196552) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298518) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(0.028129015) q[0];
rz(-1.4460538) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(0.45949724) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4079553) q[0];
sx q[0];
rz(-1.6064294) q[0];
sx q[0];
rz(0.028546988) q[0];
rz(-pi) q[1];
rz(2.8064583) q[2];
sx q[2];
rz(-2.2907933) q[2];
sx q[2];
rz(-1.1283666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34950911) q[1];
sx q[1];
rz(-2.4448423) q[1];
sx q[1];
rz(-1.7972444) q[1];
rz(1.2351843) q[3];
sx q[3];
rz(-0.40235717) q[3];
sx q[3];
rz(1.622792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10099899) q[2];
sx q[2];
rz(-1.3446292) q[2];
sx q[2];
rz(2.8864268) q[2];
rz(2.1549639) q[3];
sx q[3];
rz(-0.41472236) q[3];
sx q[3];
rz(-1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3808909) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(-1.799452) q[0];
rz(0.0019207151) q[1];
sx q[1];
rz(-1.0425967) q[1];
sx q[1];
rz(1.049918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.130803) q[0];
sx q[0];
rz(-1.783004) q[0];
sx q[0];
rz(-0.15168587) q[0];
x q[1];
rz(-1.7467612) q[2];
sx q[2];
rz(-1.5696042) q[2];
sx q[2];
rz(-0.81874412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.695076) q[1];
sx q[1];
rz(-0.43276603) q[1];
sx q[1];
rz(0.54181918) q[1];
x q[2];
rz(1.330659) q[3];
sx q[3];
rz(-0.74455816) q[3];
sx q[3];
rz(-1.9743686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2125825) q[2];
sx q[2];
rz(-1.9276103) q[2];
sx q[2];
rz(-2.6455961) q[2];
rz(-1.4604733) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7620508) q[0];
sx q[0];
rz(-1.1045781) q[0];
sx q[0];
rz(-1.7212575) q[0];
rz(-0.63915359) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(0.47472246) q[2];
sx q[2];
rz(-1.793515) q[2];
sx q[2];
rz(-1.1303177) q[2];
rz(0.81845508) q[3];
sx q[3];
rz(-1.0651121) q[3];
sx q[3];
rz(-2.4218694) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
