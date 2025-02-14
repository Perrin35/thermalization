OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.017555822) q[0];
sx q[0];
rz(3.2942009) q[0];
sx q[0];
rz(9.6080544) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(-1.1396989) q[1];
sx q[1];
rz(-0.5782063) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59435049) q[0];
sx q[0];
rz(-0.092012398) q[0];
sx q[0];
rz(0.1081744) q[0];
x q[1];
rz(1.6553641) q[2];
sx q[2];
rz(-1.5377759) q[2];
sx q[2];
rz(1.7109035) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44408724) q[1];
sx q[1];
rz(-1.1050535) q[1];
sx q[1];
rz(-2.5986141) q[1];
rz(-pi) q[2];
rz(-1.1050866) q[3];
sx q[3];
rz(-1.5381201) q[3];
sx q[3];
rz(1.0099883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1766498) q[2];
sx q[2];
rz(-2.2735333) q[2];
sx q[2];
rz(3.0435666) q[2];
rz(2.3227504) q[3];
sx q[3];
rz(-0.10418532) q[3];
sx q[3];
rz(0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0007415) q[0];
sx q[0];
rz(-0.31743693) q[0];
sx q[0];
rz(2.4733518) q[0];
rz(2.5375598) q[1];
sx q[1];
rz(-0.030345358) q[1];
sx q[1];
rz(2.277453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0111496) q[0];
sx q[0];
rz(-2.1571223) q[0];
sx q[0];
rz(-1.4545813) q[0];
rz(2.2496944) q[2];
sx q[2];
rz(-1.9111655) q[2];
sx q[2];
rz(2.3838885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64795463) q[1];
sx q[1];
rz(-1.3820547) q[1];
sx q[1];
rz(-0.24540875) q[1];
rz(0.20012466) q[3];
sx q[3];
rz(-1.311655) q[3];
sx q[3];
rz(-2.5166496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.981367) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(-2.3571864) q[2];
rz(1.9085599) q[3];
sx q[3];
rz(-0.27850702) q[3];
sx q[3];
rz(-1.2109141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3311555) q[0];
sx q[0];
rz(-1.5758608) q[0];
sx q[0];
rz(2.1203777) q[0];
rz(1.5039697) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(-0.57600299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.841678) q[0];
sx q[0];
rz(-1.9865216) q[0];
sx q[0];
rz(-1.6982416) q[0];
x q[1];
rz(0.25635135) q[2];
sx q[2];
rz(-1.1927989) q[2];
sx q[2];
rz(2.5963654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4788379) q[1];
sx q[1];
rz(-2.38677) q[1];
sx q[1];
rz(-0.50601064) q[1];
rz(-pi) q[2];
rz(-0.27184527) q[3];
sx q[3];
rz(-1.3619641) q[3];
sx q[3];
rz(-2.6110716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3699428) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(-0.096262781) q[2];
rz(-2.7104968) q[3];
sx q[3];
rz(-0.96612203) q[3];
sx q[3];
rz(2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189321) q[0];
sx q[0];
rz(-0.12598251) q[0];
sx q[0];
rz(0.2952964) q[0];
rz(2.2604306) q[1];
sx q[1];
rz(-2.9144139) q[1];
sx q[1];
rz(-0.49240246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2797564) q[0];
sx q[0];
rz(-1.4040213) q[0];
sx q[0];
rz(1.357973) q[0];
rz(-pi) q[1];
rz(-1.6642163) q[2];
sx q[2];
rz(-1.5216516) q[2];
sx q[2];
rz(1.3256595) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4223692) q[1];
sx q[1];
rz(-1.2417267) q[1];
sx q[1];
rz(-0.13767875) q[1];
rz(-0.6564859) q[3];
sx q[3];
rz(-0.53561775) q[3];
sx q[3];
rz(1.519062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.015335036) q[2];
sx q[2];
rz(-0.33053645) q[2];
sx q[2];
rz(2.7988561) q[2];
rz(0.98894173) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(-0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88687503) q[0];
sx q[0];
rz(-0.29061341) q[0];
sx q[0];
rz(1.8732204) q[0];
rz(-0.36240029) q[1];
sx q[1];
rz(-2.2762894) q[1];
sx q[1];
rz(-1.5579582) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0104631) q[0];
sx q[0];
rz(-2.4005425) q[0];
sx q[0];
rz(0.69990943) q[0];
rz(0.18581219) q[2];
sx q[2];
rz(-1.8430329) q[2];
sx q[2];
rz(0.21493064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4526738) q[1];
sx q[1];
rz(-0.59725475) q[1];
sx q[1];
rz(-3.1147472) q[1];
x q[2];
rz(2.5592119) q[3];
sx q[3];
rz(-0.60027585) q[3];
sx q[3];
rz(1.511661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6543988) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-2.2844592) q[2];
rz(-1.0355787) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(0.63151675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71984464) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(-0.33472043) q[0];
rz(-0.98182976) q[1];
sx q[1];
rz(-1.5515168) q[1];
sx q[1];
rz(-2.2614711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3405238) q[0];
sx q[0];
rz(-2.8084917) q[0];
sx q[0];
rz(2.0621501) q[0];
rz(-pi) q[1];
rz(-0.99079972) q[2];
sx q[2];
rz(-0.31337619) q[2];
sx q[2];
rz(0.51381451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33226686) q[1];
sx q[1];
rz(-2.6490445) q[1];
sx q[1];
rz(1.2713856) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2273034) q[3];
sx q[3];
rz(-0.56826545) q[3];
sx q[3];
rz(1.9075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6191787) q[2];
sx q[2];
rz(-2.652707) q[2];
sx q[2];
rz(-1.6467113) q[2];
rz(1.304168) q[3];
sx q[3];
rz(-0.84916484) q[3];
sx q[3];
rz(2.9554534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3592767) q[0];
sx q[0];
rz(-2.9501811) q[0];
sx q[0];
rz(3.0707448) q[0];
rz(0.62295667) q[1];
sx q[1];
rz(-0.40095913) q[1];
sx q[1];
rz(2.7080022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569048) q[0];
sx q[0];
rz(-1.5400209) q[0];
sx q[0];
rz(2.8367443) q[0];
rz(-pi) q[1];
x q[1];
rz(2.95288) q[2];
sx q[2];
rz(-0.29653877) q[2];
sx q[2];
rz(1.6558356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63804352) q[1];
sx q[1];
rz(-1.2473599) q[1];
sx q[1];
rz(1.6721318) q[1];
rz(2.6979792) q[3];
sx q[3];
rz(-1.6742236) q[3];
sx q[3];
rz(-0.87956968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5101461) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(0.637429) q[2];
rz(-2.0006477) q[3];
sx q[3];
rz(-2.7009522) q[3];
sx q[3];
rz(2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7618074) q[0];
sx q[0];
rz(-0.26476911) q[0];
sx q[0];
rz(-2.8763212) q[0];
rz(3.1130262) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(-2.2144337) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521479) q[0];
sx q[0];
rz(-1.5434221) q[0];
sx q[0];
rz(0.22449799) q[0];
rz(-0.36225944) q[2];
sx q[2];
rz(-2.7484244) q[2];
sx q[2];
rz(2.071381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88419848) q[1];
sx q[1];
rz(-0.78620877) q[1];
sx q[1];
rz(1.4255852) q[1];
rz(-pi) q[2];
rz(-2.718247) q[3];
sx q[3];
rz(-1.6739427) q[3];
sx q[3];
rz(2.1684017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11904968) q[2];
sx q[2];
rz(-1.0485342) q[2];
sx q[2];
rz(-2.0691464) q[2];
rz(0.33506814) q[3];
sx q[3];
rz(-3.0266893) q[3];
sx q[3];
rz(-0.2255628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91113126) q[0];
sx q[0];
rz(-1.1536396) q[0];
sx q[0];
rz(-0.78900868) q[0];
rz(2.543653) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(-2.423563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563468) q[0];
sx q[0];
rz(-1.5296773) q[0];
sx q[0];
rz(1.4151735) q[0];
x q[1];
rz(-1.6512376) q[2];
sx q[2];
rz(-3.0249964) q[2];
sx q[2];
rz(1.2355905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8666229) q[1];
sx q[1];
rz(-0.94511813) q[1];
sx q[1];
rz(0.83557202) q[1];
x q[2];
rz(-1.6634462) q[3];
sx q[3];
rz(-1.3934373) q[3];
sx q[3];
rz(-2.5611344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.614552) q[2];
sx q[2];
rz(-0.25200945) q[2];
sx q[2];
rz(1.7238114) q[2];
rz(-2.0333911) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(-0.38780701) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.336816) q[0];
sx q[0];
rz(-2.5122061) q[0];
sx q[0];
rz(0.25548536) q[0];
rz(-0.77087036) q[1];
sx q[1];
rz(-1.5893385) q[1];
sx q[1];
rz(1.593387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436018) q[0];
sx q[0];
rz(-1.457347) q[0];
sx q[0];
rz(-1.1313637) q[0];
rz(-pi) q[1];
rz(-0.86779197) q[2];
sx q[2];
rz(-1.7475339) q[2];
sx q[2];
rz(2.8121626) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0615427) q[1];
sx q[1];
rz(-1.1876186) q[1];
sx q[1];
rz(1.2754759) q[1];
x q[2];
rz(-2.6029165) q[3];
sx q[3];
rz(-0.65924257) q[3];
sx q[3];
rz(1.074312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-2.3994583) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(0.54002386) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(-2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.222432) q[0];
sx q[0];
rz(-1.5615015) q[0];
sx q[0];
rz(-1.1175565) q[0];
rz(0.22668214) q[1];
sx q[1];
rz(-3.0381028) q[1];
sx q[1];
rz(1.4729952) q[1];
rz(2.1508455) q[2];
sx q[2];
rz(-1.8875835) q[2];
sx q[2];
rz(0.86474309) q[2];
rz(-1.8899824) q[3];
sx q[3];
rz(-0.87379119) q[3];
sx q[3];
rz(-2.9004267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
