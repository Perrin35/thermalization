OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3992231) q[0];
sx q[0];
rz(-2.7288781) q[0];
sx q[0];
rz(-0.8362008) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(1.0083415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158186) q[0];
sx q[0];
rz(-1.4918707) q[0];
sx q[0];
rz(-1.1962587) q[0];
rz(-pi) q[1];
rz(0.54097367) q[2];
sx q[2];
rz(-0.43593513) q[2];
sx q[2];
rz(-1.774314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.28532449) q[1];
sx q[1];
rz(-0.50438499) q[1];
sx q[1];
rz(2.8125202) q[1];
rz(1.3775695) q[3];
sx q[3];
rz(-1.591103) q[3];
sx q[3];
rz(-1.3042252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3081554) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(2.4885139) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(1.0623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(-1.1603181) q[1];
sx q[1];
rz(-2.7313045) q[1];
sx q[1];
rz(0.62058273) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4575093) q[0];
sx q[0];
rz(-0.39827049) q[0];
sx q[0];
rz(-1.0578367) q[0];
rz(-pi) q[1];
rz(-2.0257904) q[2];
sx q[2];
rz(-2.4206851) q[2];
sx q[2];
rz(-1.7146595) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0936135) q[1];
sx q[1];
rz(-1.9318214) q[1];
sx q[1];
rz(0.71012902) q[1];
rz(-2.1866131) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(-1.3835088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6985942) q[2];
sx q[2];
rz(-1.5195165) q[2];
sx q[2];
rz(-2.911574) q[2];
rz(-1.5864141) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(-2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25552937) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(-2.0607167) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(0.08531514) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9935308) q[0];
sx q[0];
rz(-0.9328273) q[0];
sx q[0];
rz(-0.83420269) q[0];
rz(-pi) q[1];
rz(1.7833461) q[2];
sx q[2];
rz(-1.0091127) q[2];
sx q[2];
rz(-1.4457653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18135611) q[1];
sx q[1];
rz(-1.4048368) q[1];
sx q[1];
rz(-1.7039429) q[1];
rz(0.063488678) q[3];
sx q[3];
rz(-2.5194019) q[3];
sx q[3];
rz(-2.6486462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9366511) q[2];
sx q[2];
rz(-0.0095657883) q[2];
sx q[2];
rz(2.9452475) q[2];
rz(-2.8593072) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(-2.096368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008217) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-0.38544449) q[0];
rz(1.7281945) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(-0.24994303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1544558) q[0];
sx q[0];
rz(-0.76351316) q[0];
sx q[0];
rz(-2.187243) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0536575) q[2];
sx q[2];
rz(-1.5107949) q[2];
sx q[2];
rz(-1.1052002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.476772) q[1];
sx q[1];
rz(-0.52174924) q[1];
sx q[1];
rz(-2.797965) q[1];
rz(-pi) q[2];
rz(1.7697261) q[3];
sx q[3];
rz(-2.3545111) q[3];
sx q[3];
rz(0.50686344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6803153) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(-2.9052022) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(-2.6564823) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(2.3050883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77959767) q[0];
sx q[0];
rz(-1.2149724) q[0];
sx q[0];
rz(-1.9076288) q[0];
x q[1];
rz(-2.2274687) q[2];
sx q[2];
rz(-2.6466922) q[2];
sx q[2];
rz(2.6656541) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5768535) q[1];
sx q[1];
rz(-0.67138636) q[1];
sx q[1];
rz(1.6480084) q[1];
rz(-2.9936166) q[3];
sx q[3];
rz(-2.736909) q[3];
sx q[3];
rz(-1.2860822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89852077) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(-2.3720429) q[2];
rz(0.92929333) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(2.9088959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-2.6564044) q[0];
sx q[0];
rz(2.3727681) q[0];
rz(1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(-2.2519055) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2425491) q[0];
sx q[0];
rz(-2.7022323) q[0];
sx q[0];
rz(1.0602632) q[0];
x q[1];
rz(1.5264411) q[2];
sx q[2];
rz(-1.5694478) q[2];
sx q[2];
rz(1.1464455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1353768) q[1];
sx q[1];
rz(-1.6172503) q[1];
sx q[1];
rz(-1.4245288) q[1];
rz(-1.6857288) q[3];
sx q[3];
rz(-2.0191666) q[3];
sx q[3];
rz(-1.2796206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49030226) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5865405) q[2];
rz(1.8709315) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(3.0095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1726058) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(-0.78713083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2427776) q[0];
sx q[0];
rz(-1.4393596) q[0];
sx q[0];
rz(-1.7064894) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93573715) q[2];
sx q[2];
rz(-0.79197403) q[2];
sx q[2];
rz(0.21312772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63185502) q[1];
sx q[1];
rz(-0.69586772) q[1];
sx q[1];
rz(0.79001714) q[1];
rz(0.58148099) q[3];
sx q[3];
rz(-0.6430917) q[3];
sx q[3];
rz(3.0742925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0234915) q[2];
sx q[2];
rz(-0.96994895) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(-1.7484885) q[0];
rz(-0.45685592) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(-2.0410062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7129214) q[0];
sx q[0];
rz(-1.572519) q[0];
sx q[0];
rz(-0.07600204) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4462561) q[2];
sx q[2];
rz(-2.0954663) q[2];
sx q[2];
rz(2.6161043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3708027) q[1];
sx q[1];
rz(-2.2375771) q[1];
sx q[1];
rz(-0.043170269) q[1];
rz(-pi) q[2];
rz(2.8131717) q[3];
sx q[3];
rz(-0.41116086) q[3];
sx q[3];
rz(2.823165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0988203) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5808606) q[2];
rz(2.3780195) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(-1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(2.8501046) q[0];
rz(-1.905722) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(-3.018697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50704573) q[0];
sx q[0];
rz(-1.5193918) q[0];
sx q[0];
rz(2.794751) q[0];
x q[1];
rz(-0.22237402) q[2];
sx q[2];
rz(-1.2952598) q[2];
sx q[2];
rz(-1.0433407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.025561573) q[1];
sx q[1];
rz(-1.4605547) q[1];
sx q[1];
rz(1.362839) q[1];
x q[2];
rz(1.6121665) q[3];
sx q[3];
rz(-2.9968516) q[3];
sx q[3];
rz(-2.5998757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5128532) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(-2.0780308) q[2];
rz(-2.9641446) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(2.1232846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43750957) q[0];
sx q[0];
rz(-0.45314416) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(2.5408632) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(-2.5240555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7970456) q[0];
sx q[0];
rz(-0.24976845) q[0];
sx q[0];
rz(2.0956371) q[0];
x q[1];
rz(-2.3138758) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(0.043443505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3869218) q[1];
sx q[1];
rz(-0.53334177) q[1];
sx q[1];
rz(-1.6120595) q[1];
x q[2];
rz(1.2631868) q[3];
sx q[3];
rz(-1.0196536) q[3];
sx q[3];
rz(1.2778417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(2.0591002) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33901535) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(0.045724178) q[1];
sx q[1];
rz(-0.8174236) q[1];
sx q[1];
rz(1.5541706) q[1];
rz(2.1004213) q[2];
sx q[2];
rz(-1.1239792) q[2];
sx q[2];
rz(-1.299528) q[2];
rz(1.4010728) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
