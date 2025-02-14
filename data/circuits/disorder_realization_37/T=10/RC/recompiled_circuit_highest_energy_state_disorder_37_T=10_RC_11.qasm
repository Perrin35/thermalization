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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(-0.65161172) q[1];
sx q[1];
rz(-1.7348644) q[1];
sx q[1];
rz(-2.9098517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.988294) q[0];
sx q[0];
rz(-0.35044119) q[0];
sx q[0];
rz(-0.20163433) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.628078) q[2];
sx q[2];
rz(-1.0250895) q[2];
sx q[2];
rz(1.6451665) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5197088) q[1];
sx q[1];
rz(-0.95920282) q[1];
sx q[1];
rz(1.1033464) q[1];
x q[2];
rz(-1.6123338) q[3];
sx q[3];
rz(-1.9088863) q[3];
sx q[3];
rz(1.5134144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0509402) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(0.88979641) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-0.92581785) q[0];
rz(2.0766808) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(-0.085478641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8259698) q[0];
sx q[0];
rz(-1.96061) q[0];
sx q[0];
rz(-1.5743544) q[0];
x q[1];
rz(-0.95909848) q[2];
sx q[2];
rz(-1.2818205) q[2];
sx q[2];
rz(-1.1346357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2457471) q[1];
sx q[1];
rz(-1.3407882) q[1];
sx q[1];
rz(1.0279015) q[1];
rz(-pi) q[2];
rz(-2.0382337) q[3];
sx q[3];
rz(-0.93884838) q[3];
sx q[3];
rz(-1.2480717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3257137) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(-0.386664) q[2];
rz(1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81095186) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(0.49705848) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(-2.846948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9884856) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(1.6436623) q[0];
rz(1.398009) q[2];
sx q[2];
rz(-0.9388939) q[2];
sx q[2];
rz(-1.2829219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68874796) q[1];
sx q[1];
rz(-1.1383245) q[1];
sx q[1];
rz(2.2562863) q[1];
rz(-pi) q[2];
rz(-1.7792652) q[3];
sx q[3];
rz(-2.5304171) q[3];
sx q[3];
rz(1.1217505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7872539) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5798689) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314986) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.6353105) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(-0.062072676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6784883) q[0];
sx q[0];
rz(-2.0651428) q[0];
sx q[0];
rz(2.2114179) q[0];
x q[1];
rz(1.9628635) q[2];
sx q[2];
rz(-1.829756) q[2];
sx q[2];
rz(2.1990537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.049438795) q[1];
sx q[1];
rz(-0.88969943) q[1];
sx q[1];
rz(-1.9309631) q[1];
rz(-pi) q[2];
rz(2.7113879) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(-1.9752432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0486003) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(-1.3108866) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(-2.766975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(2.2485961) q[0];
rz(2.112174) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(0.095887862) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3678148) q[0];
sx q[0];
rz(-2.4565963) q[0];
sx q[0];
rz(0.8789341) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0711543) q[2];
sx q[2];
rz(-1.4820921) q[2];
sx q[2];
rz(-2.8005637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7815302) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(2.5468154) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5378102) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(-2.3100738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-0.42547697) q[2];
rz(-0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285889) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(1.0711627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334918) q[0];
sx q[0];
rz(-2.0530434) q[0];
sx q[0];
rz(1.2511176) q[0];
rz(1.2920678) q[2];
sx q[2];
rz(-1.0459002) q[2];
sx q[2];
rz(-0.43290813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6630604) q[1];
sx q[1];
rz(-2.2473865) q[1];
sx q[1];
rz(-2.8512136) q[1];
x q[2];
rz(2.8240194) q[3];
sx q[3];
rz(-1.4151207) q[3];
sx q[3];
rz(-0.50567108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57227197) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(-0.039483698) q[2];
rz(3.0651921) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.9901265) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(-1.8340402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261953) q[0];
sx q[0];
rz(-1.2457677) q[0];
sx q[0];
rz(0.23120489) q[0];
rz(-pi) q[1];
rz(1.4979532) q[2];
sx q[2];
rz(-1.9460287) q[2];
sx q[2];
rz(-2.2323687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7286018) q[1];
sx q[1];
rz(-0.54122335) q[1];
sx q[1];
rz(1.1515929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7634972) q[3];
sx q[3];
rz(-2.5268838) q[3];
sx q[3];
rz(2.4845882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82039708) q[2];
sx q[2];
rz(-0.68009192) q[2];
sx q[2];
rz(0.55366984) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(-1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(-0.30690646) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(-1.2409522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4590603) q[0];
sx q[0];
rz(-2.0143565) q[0];
sx q[0];
rz(0.17291594) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1419134) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(2.6311324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1859295) q[1];
sx q[1];
rz(-0.36500442) q[1];
sx q[1];
rz(-2.5899653) q[1];
x q[2];
rz(2.5316174) q[3];
sx q[3];
rz(-1.2362567) q[3];
sx q[3];
rz(1.1776678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79157311) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.2903068) q[2];
rz(2.4604515) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.5017989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27555585) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(1.1736897) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(-0.079708286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3143334) q[0];
sx q[0];
rz(-1.9406576) q[0];
sx q[0];
rz(2.8989559) q[0];
rz(-pi) q[1];
rz(2.9201512) q[2];
sx q[2];
rz(-1.6802854) q[2];
sx q[2];
rz(-1.2564645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6882872) q[1];
sx q[1];
rz(-0.88677471) q[1];
sx q[1];
rz(-0.75835336) q[1];
rz(-pi) q[2];
x q[2];
rz(2.856273) q[3];
sx q[3];
rz(-2.3148708) q[3];
sx q[3];
rz(-0.43275012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6012663) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(1.3339174) q[2];
rz(-1.5909083) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(-2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48225668) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.2905066) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528906) q[0];
sx q[0];
rz(-1.9963309) q[0];
sx q[0];
rz(-2.8136926) q[0];
rz(-pi) q[1];
rz(0.89628156) q[2];
sx q[2];
rz(-2.7283629) q[2];
sx q[2];
rz(0.21597029) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86856213) q[1];
sx q[1];
rz(-2.2924463) q[1];
sx q[1];
rz(2.9356586) q[1];
rz(1.4560111) q[3];
sx q[3];
rz(-1.2437395) q[3];
sx q[3];
rz(2.9714874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(1.8168137) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751547) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(-1.4229763) q[1];
sx q[1];
rz(-1.9444793) q[1];
sx q[1];
rz(0.3107298) q[1];
rz(1.995261) q[2];
sx q[2];
rz(-0.88256114) q[2];
sx q[2];
rz(0.12086856) q[2];
rz(-2.041009) q[3];
sx q[3];
rz(-0.96405021) q[3];
sx q[3];
rz(-1.8586803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
