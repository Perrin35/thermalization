OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.032441) q[0];
sx q[0];
rz(-1.1530131) q[0];
sx q[0];
rz(3.0670526) q[0];
rz(1.6115161) q[1];
sx q[1];
rz(3.2009701) q[1];
sx q[1];
rz(8.9006807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7439047) q[0];
sx q[0];
rz(-0.70262229) q[0];
sx q[0];
rz(0.98708679) q[0];
rz(2.5251389) q[2];
sx q[2];
rz(-0.40531593) q[2];
sx q[2];
rz(-2.7242994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9896696) q[1];
sx q[1];
rz(-2.1300585) q[1];
sx q[1];
rz(-0.86821809) q[1];
rz(-1.662781) q[3];
sx q[3];
rz(-1.5301203) q[3];
sx q[3];
rz(-2.3856861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3268299) q[2];
sx q[2];
rz(-0.49850285) q[2];
sx q[2];
rz(-0.89547431) q[2];
rz(0.26252663) q[3];
sx q[3];
rz(-2.7956876) q[3];
sx q[3];
rz(3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9240016) q[0];
sx q[0];
rz(-1.1955248) q[0];
sx q[0];
rz(-0.1161282) q[0];
rz(0.63236347) q[1];
sx q[1];
rz(-2.875681) q[1];
sx q[1];
rz(-1.6557453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350548) q[0];
sx q[0];
rz(-3.0710906) q[0];
sx q[0];
rz(-2.8126731) q[0];
rz(0.2367649) q[2];
sx q[2];
rz(-2.1000266) q[2];
sx q[2];
rz(2.6108612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8101059) q[1];
sx q[1];
rz(-1.4432194) q[1];
sx q[1];
rz(1.4639616) q[1];
rz(-2.2712098) q[3];
sx q[3];
rz(-1.4705949) q[3];
sx q[3];
rz(-2.9068275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8302725) q[2];
sx q[2];
rz(-0.98036426) q[2];
sx q[2];
rz(-0.92973989) q[2];
rz(-2.7021507) q[3];
sx q[3];
rz(-1.6012871) q[3];
sx q[3];
rz(-0.75696993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0759401) q[0];
sx q[0];
rz(-0.75242281) q[0];
sx q[0];
rz(0.886985) q[0];
rz(-1.8807962) q[1];
sx q[1];
rz(-2.0340684) q[1];
sx q[1];
rz(-1.4874123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9943574) q[0];
sx q[0];
rz(-1.571313) q[0];
sx q[0];
rz(-1.4532671) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4753129) q[2];
sx q[2];
rz(-1.9621358) q[2];
sx q[2];
rz(-0.0454768) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9076242) q[1];
sx q[1];
rz(-0.32999906) q[1];
sx q[1];
rz(-0.054363175) q[1];
rz(-pi) q[2];
rz(2.5599407) q[3];
sx q[3];
rz(-1.8466788) q[3];
sx q[3];
rz(-1.1946071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6416574) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(0.43528834) q[2];
rz(-0.38517243) q[3];
sx q[3];
rz(-1.2110854) q[3];
sx q[3];
rz(-0.96863532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73548549) q[0];
sx q[0];
rz(-0.084097363) q[0];
sx q[0];
rz(-0.054585833) q[0];
rz(1.5054043) q[1];
sx q[1];
rz(-1.9278229) q[1];
sx q[1];
rz(2.7240567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62359257) q[0];
sx q[0];
rz(-1.5125649) q[0];
sx q[0];
rz(0.43727711) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63234858) q[2];
sx q[2];
rz(-1.8397619) q[2];
sx q[2];
rz(-0.40695813) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3247973) q[1];
sx q[1];
rz(-0.15843219) q[1];
sx q[1];
rz(0.043912391) q[1];
rz(-pi) q[2];
rz(-1.8330386) q[3];
sx q[3];
rz(-1.635437) q[3];
sx q[3];
rz(1.3418947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85228449) q[2];
sx q[2];
rz(-0.70312971) q[2];
sx q[2];
rz(-0.0011778041) q[2];
rz(-1.3834472) q[3];
sx q[3];
rz(-0.26418424) q[3];
sx q[3];
rz(-2.0980825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99172878) q[0];
sx q[0];
rz(-2.9809451) q[0];
sx q[0];
rz(2.7657261) q[0];
rz(2.1916892) q[1];
sx q[1];
rz(-1.6641195) q[1];
sx q[1];
rz(-0.72521597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26869088) q[0];
sx q[0];
rz(-0.16624545) q[0];
sx q[0];
rz(-2.3670271) q[0];
x q[1];
rz(2.2711168) q[2];
sx q[2];
rz(-1.0511802) q[2];
sx q[2];
rz(-0.82110559) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2754073) q[1];
sx q[1];
rz(-1.5650041) q[1];
sx q[1];
rz(-0.097472982) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3106724) q[3];
sx q[3];
rz(-0.88444607) q[3];
sx q[3];
rz(-0.42410606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3788562) q[2];
sx q[2];
rz(-1.4521658) q[2];
sx q[2];
rz(2.6689996) q[2];
rz(3.0224814) q[3];
sx q[3];
rz(-2.5178858) q[3];
sx q[3];
rz(-0.81010336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3177719) q[0];
sx q[0];
rz(-2.9404984) q[0];
sx q[0];
rz(0.066545181) q[0];
rz(-2.9464974) q[1];
sx q[1];
rz(-1.0761484) q[1];
sx q[1];
rz(1.9871575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46259743) q[0];
sx q[0];
rz(-0.53946153) q[0];
sx q[0];
rz(-2.0360002) q[0];
x q[1];
rz(-1.4835266) q[2];
sx q[2];
rz(-1.1157398) q[2];
sx q[2];
rz(-0.94489646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2278359) q[1];
sx q[1];
rz(-1.933473) q[1];
sx q[1];
rz(-3.0272724) q[1];
x q[2];
rz(1.5905753) q[3];
sx q[3];
rz(-0.58073509) q[3];
sx q[3];
rz(1.6466717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92516148) q[2];
sx q[2];
rz(-1.5387646) q[2];
sx q[2];
rz(0.72539854) q[2];
rz(1.7321436) q[3];
sx q[3];
rz(-1.9603445) q[3];
sx q[3];
rz(0.60666549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055939097) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(-0.90580171) q[0];
rz(0.55465758) q[1];
sx q[1];
rz(-0.83811086) q[1];
sx q[1];
rz(-2.99756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54948037) q[0];
sx q[0];
rz(-1.0735816) q[0];
sx q[0];
rz(-0.29659941) q[0];
x q[1];
rz(0.085113581) q[2];
sx q[2];
rz(-2.0817698) q[2];
sx q[2];
rz(2.2764595) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0195469) q[1];
sx q[1];
rz(-2.658469) q[1];
sx q[1];
rz(-2.2058317) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1062076) q[3];
sx q[3];
rz(-0.28942063) q[3];
sx q[3];
rz(1.1076224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4732699) q[2];
sx q[2];
rz(-2.747135) q[2];
sx q[2];
rz(-0.99297601) q[2];
rz(0.18443491) q[3];
sx q[3];
rz(-0.95576972) q[3];
sx q[3];
rz(-0.9930281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011878012) q[0];
sx q[0];
rz(-0.60420245) q[0];
sx q[0];
rz(-0.39464828) q[0];
rz(-0.53798419) q[1];
sx q[1];
rz(-2.4514276) q[1];
sx q[1];
rz(-1.1916377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906843) q[0];
sx q[0];
rz(-2.5328851) q[0];
sx q[0];
rz(2.6842791) q[0];
rz(-pi) q[1];
rz(-0.46311867) q[2];
sx q[2];
rz(-0.41234932) q[2];
sx q[2];
rz(-0.3316484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.887693) q[1];
sx q[1];
rz(-2.116345) q[1];
sx q[1];
rz(2.6274526) q[1];
x q[2];
rz(1.1147145) q[3];
sx q[3];
rz(-0.9852162) q[3];
sx q[3];
rz(2.8600313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3523606) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(2.1996876) q[2];
rz(0.76147979) q[3];
sx q[3];
rz(-1.1250291) q[3];
sx q[3];
rz(3.1010845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86904675) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(-2.0696562) q[0];
rz(0.6140703) q[1];
sx q[1];
rz(-0.89659381) q[1];
sx q[1];
rz(0.44642064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9337685) q[0];
sx q[0];
rz(-2.38777) q[0];
sx q[0];
rz(-1.1541744) q[0];
rz(-pi) q[1];
rz(-2.1062966) q[2];
sx q[2];
rz(-0.95343325) q[2];
sx q[2];
rz(1.8445002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.804396) q[1];
sx q[1];
rz(-2.2257518) q[1];
sx q[1];
rz(1.7095997) q[1];
rz(-pi) q[2];
rz(-0.50531252) q[3];
sx q[3];
rz(-1.4243505) q[3];
sx q[3];
rz(0.088602655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3093695) q[2];
sx q[2];
rz(-2.1840405) q[2];
sx q[2];
rz(1.2584125) q[2];
rz(-1.2331412) q[3];
sx q[3];
rz(-0.62658739) q[3];
sx q[3];
rz(-1.3174177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719139) q[0];
sx q[0];
rz(-0.8154251) q[0];
sx q[0];
rz(1.7183787) q[0];
rz(-2.8990959) q[1];
sx q[1];
rz(-0.47912326) q[1];
sx q[1];
rz(1.2538145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72287382) q[0];
sx q[0];
rz(-0.44298816) q[0];
sx q[0];
rz(1.1305869) q[0];
x q[1];
rz(-1.5635339) q[2];
sx q[2];
rz(-1.6374131) q[2];
sx q[2];
rz(2.3965701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1088646) q[1];
sx q[1];
rz(-1.0582339) q[1];
sx q[1];
rz(-2.7973919) q[1];
rz(-0.21988545) q[3];
sx q[3];
rz(-1.2078834) q[3];
sx q[3];
rz(-2.9232034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32580507) q[2];
sx q[2];
rz(-0.82745224) q[2];
sx q[2];
rz(-3.117756) q[2];
rz(1.2787974) q[3];
sx q[3];
rz(-0.33134225) q[3];
sx q[3];
rz(2.3648025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9780289) q[0];
sx q[0];
rz(-1.4373056) q[0];
sx q[0];
rz(1.9594255) q[0];
rz(0.72721807) q[1];
sx q[1];
rz(-1.4677508) q[1];
sx q[1];
rz(-0.8263091) q[1];
rz(0.019451774) q[2];
sx q[2];
rz(-1.7218334) q[2];
sx q[2];
rz(-1.3109372) q[2];
rz(-1.1990697) q[3];
sx q[3];
rz(-3.0109497) q[3];
sx q[3];
rz(2.0554832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
