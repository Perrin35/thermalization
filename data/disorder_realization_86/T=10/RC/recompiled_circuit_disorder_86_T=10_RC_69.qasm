OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1887814) q[0];
sx q[0];
rz(-0.45438284) q[0];
sx q[0];
rz(-0.36034583) q[0];
rz(-pi) q[1];
rz(-1.7035159) q[2];
sx q[2];
rz(-1.3157529) q[2];
sx q[2];
rz(-2.1357352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23956242) q[1];
sx q[1];
rz(-2.4671577) q[1];
sx q[1];
rz(-2.2873138) q[1];
x q[2];
rz(-3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(-1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(0.96347934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.11461) q[0];
sx q[0];
rz(-0.61404213) q[0];
sx q[0];
rz(1.5668037) q[0];
rz(-2.8112667) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(2.3842173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0429093) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(0.98867464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9217334) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(0.79743687) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(-0.55999666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277498) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(-1.9613128) q[0];
rz(2.838344) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(0.081239935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3126038) q[1];
sx q[1];
rz(-1.7692411) q[1];
sx q[1];
rz(2.7752084) q[1];
rz(2.3476944) q[3];
sx q[3];
rz(-2.5069935) q[3];
sx q[3];
rz(1.3018228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(-0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(-1.8428615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7553058) q[0];
sx q[0];
rz(-2.5884429) q[0];
sx q[0];
rz(-2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5694593) q[2];
sx q[2];
rz(-0.7664116) q[2];
sx q[2];
rz(-3.1198451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97949308) q[1];
sx q[1];
rz(-2.5433308) q[1];
sx q[1];
rz(1.9059327) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2508568) q[3];
sx q[3];
rz(-0.3443998) q[3];
sx q[3];
rz(0.26564769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5359042) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(-2.9105913) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(-2.8447661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2192229) q[0];
sx q[0];
rz(-2.5693359) q[0];
sx q[0];
rz(-3.0692283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9752521) q[2];
sx q[2];
rz(-0.74976774) q[2];
sx q[2];
rz(1.9094085) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27767147) q[1];
sx q[1];
rz(-1.1754981) q[1];
sx q[1];
rz(-1.0635832) q[1];
rz(-pi) q[2];
rz(-2.4279168) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(-1.0825368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(2.3902067) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-2.5352535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6316815) q[0];
sx q[0];
rz(-1.5874377) q[0];
sx q[0];
rz(-3.095754) q[0];
x q[1];
rz(2.2491127) q[2];
sx q[2];
rz(-1.9024897) q[2];
sx q[2];
rz(1.4075556) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(-1.238766) q[1];
x q[2];
rz(-1.7883676) q[3];
sx q[3];
rz(-1.4975582) q[3];
sx q[3];
rz(2.6069802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(2.3513444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761844) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(2.5248812) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21364084) q[2];
sx q[2];
rz(-1.5593312) q[2];
sx q[2];
rz(-1.2917047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.33752791) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(2.7359664) q[1];
x q[2];
rz(-2.4942057) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(-3.0800381) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247894) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.75288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6435218) q[0];
sx q[0];
rz(-2.5812491) q[0];
sx q[0];
rz(1.5146921) q[0];
rz(-pi) q[1];
rz(1.9954761) q[2];
sx q[2];
rz(-1.0815902) q[2];
sx q[2];
rz(-2.3236772) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(-3.0180879) q[1];
rz(2.8453313) q[3];
sx q[3];
rz(-0.98516609) q[3];
sx q[3];
rz(-1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(-0.88225538) q[2];
rz(-1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(2.5833599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840069) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(0.42752479) q[0];
rz(0.054106601) q[2];
sx q[2];
rz(-1.4634842) q[2];
sx q[2];
rz(3.0215614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0135632) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(-3.0548884) q[1];
x q[2];
rz(-1.3853119) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(-0.51717796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1062746) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(-1.6465228) q[0];
rz(3.0977544) q[2];
sx q[2];
rz(-2.1423116) q[2];
sx q[2];
rz(0.29585719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2460675) q[1];
sx q[1];
rz(-1.2125891) q[1];
sx q[1];
rz(-1.9073061) q[1];
x q[2];
rz(-1.80199) q[3];
sx q[3];
rz(-1.8095067) q[3];
sx q[3];
rz(-1.064144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-3.1303828) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(2.3446463) q[3];
sx q[3];
rz(-1.7399825) q[3];
sx q[3];
rz(0.83132838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
