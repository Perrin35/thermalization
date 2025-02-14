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
rz(1.1462829) q[0];
sx q[0];
rz(-3.1132071) q[0];
sx q[0];
rz(-0.95771587) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(-1.6038916) q[1];
sx q[1];
rz(0.62827194) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.385163) q[0];
sx q[0];
rz(-1.9397501) q[0];
sx q[0];
rz(-1.4403309) q[0];
x q[1];
rz(-0.45189721) q[2];
sx q[2];
rz(-2.0263645) q[2];
sx q[2];
rz(-1.4939552) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4406887) q[1];
sx q[1];
rz(-1.2504185) q[1];
sx q[1];
rz(0.87030555) q[1];
x q[2];
rz(1.98313) q[3];
sx q[3];
rz(-1.6618007) q[3];
sx q[3];
rz(3.0418398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4665224) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(-2.7543219) q[2];
rz(1.1747423) q[3];
sx q[3];
rz(-1.4459123) q[3];
sx q[3];
rz(0.89850473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627477) q[0];
sx q[0];
rz(-1.7411106) q[0];
sx q[0];
rz(-1.8722906) q[0];
rz(1.4461888) q[1];
sx q[1];
rz(-1.2934928) q[1];
sx q[1];
rz(-0.92099014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5808153) q[0];
sx q[0];
rz(-1.2755503) q[0];
sx q[0];
rz(-0.92714374) q[0];
rz(1.3878421) q[2];
sx q[2];
rz(-0.91380807) q[2];
sx q[2];
rz(2.8859101) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16352902) q[1];
sx q[1];
rz(-2.7394751) q[1];
sx q[1];
rz(1.0671713) q[1];
rz(-2.9156906) q[3];
sx q[3];
rz(-1.5776565) q[3];
sx q[3];
rz(3.1202735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9903119) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(0.82378236) q[2];
rz(-1.524205) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(-1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9059432) q[0];
sx q[0];
rz(-0.88053572) q[0];
sx q[0];
rz(-2.549951) q[0];
rz(0.94332424) q[1];
sx q[1];
rz(-1.0572409) q[1];
sx q[1];
rz(-1.0234157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9647261) q[0];
sx q[0];
rz(-0.15892488) q[0];
sx q[0];
rz(2.311241) q[0];
x q[1];
rz(1.0092956) q[2];
sx q[2];
rz(-0.63234419) q[2];
sx q[2];
rz(0.48582669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3262965) q[1];
sx q[1];
rz(-3.0440512) q[1];
sx q[1];
rz(-1.8952712) q[1];
rz(-pi) q[2];
rz(-0.68444201) q[3];
sx q[3];
rz(-1.8307643) q[3];
sx q[3];
rz(0.79567676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8817875) q[2];
sx q[2];
rz(-1.9116348) q[2];
sx q[2];
rz(1.7477431) q[2];
rz(1.7143837) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(-0.29903308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77089906) q[0];
sx q[0];
rz(-0.40632668) q[0];
sx q[0];
rz(-2.4942177) q[0];
rz(0.24230832) q[1];
sx q[1];
rz(-1.7055885) q[1];
sx q[1];
rz(1.4586331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3998916) q[0];
sx q[0];
rz(-1.2477324) q[0];
sx q[0];
rz(1.7123651) q[0];
rz(2.5299923) q[2];
sx q[2];
rz(-1.0047874) q[2];
sx q[2];
rz(1.5274439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35136552) q[1];
sx q[1];
rz(-1.0855165) q[1];
sx q[1];
rz(-0.041260551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2599873) q[3];
sx q[3];
rz(-2.2223916) q[3];
sx q[3];
rz(2.6736265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4316537) q[2];
sx q[2];
rz(-2.2971051) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(2.8905408) q[3];
sx q[3];
rz(-2.4243088) q[3];
sx q[3];
rz(2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5421903) q[0];
sx q[0];
rz(-1.1325862) q[0];
sx q[0];
rz(-0.19790025) q[0];
rz(1.2359765) q[1];
sx q[1];
rz(-2.7695152) q[1];
sx q[1];
rz(3.1370251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3749927) q[0];
sx q[0];
rz(-2.1439634) q[0];
sx q[0];
rz(-1.8660924) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7103689) q[2];
sx q[2];
rz(-0.53832952) q[2];
sx q[2];
rz(2.0135422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98036042) q[1];
sx q[1];
rz(-1.8440108) q[1];
sx q[1];
rz(-1.7262162) q[1];
rz(-pi) q[2];
rz(2.6130535) q[3];
sx q[3];
rz(-0.40088568) q[3];
sx q[3];
rz(1.8141845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4602451) q[2];
sx q[2];
rz(-2.6474417) q[2];
sx q[2];
rz(0.34026185) q[2];
rz(1.6081238) q[3];
sx q[3];
rz(-1.7483277) q[3];
sx q[3];
rz(-1.5197915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.915864) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(-1.9541784) q[0];
rz(1.4215218) q[1];
sx q[1];
rz(-2.7233796) q[1];
sx q[1];
rz(-3.0112867) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7447259) q[0];
sx q[0];
rz(-1.6292185) q[0];
sx q[0];
rz(1.4471884) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2557507) q[2];
sx q[2];
rz(-2.3315213) q[2];
sx q[2];
rz(2.5145234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2060616) q[1];
sx q[1];
rz(-1.5714679) q[1];
sx q[1];
rz(2.9737581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8335329) q[3];
sx q[3];
rz(-1.5612226) q[3];
sx q[3];
rz(0.75754091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3336031) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(2.1017334) q[2];
rz(-0.17397675) q[3];
sx q[3];
rz(-1.1287374) q[3];
sx q[3];
rz(-2.4719293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.533605) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(2.3002891) q[0];
rz(1.3249506) q[1];
sx q[1];
rz(-2.1957896) q[1];
sx q[1];
rz(1.4680877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8155839) q[0];
sx q[0];
rz(-1.9202613) q[0];
sx q[0];
rz(-0.38527003) q[0];
rz(-3.1300342) q[2];
sx q[2];
rz(-2.0894139) q[2];
sx q[2];
rz(-0.29588884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3716814) q[1];
sx q[1];
rz(-1.103319) q[1];
sx q[1];
rz(2.5073696) q[1];
rz(-pi) q[2];
rz(-1.928059) q[3];
sx q[3];
rz(-1.284984) q[3];
sx q[3];
rz(-0.88102007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2738721) q[2];
sx q[2];
rz(-2.2404859) q[2];
sx q[2];
rz(2.4967398) q[2];
rz(-2.0598748) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(-2.4206415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69407392) q[0];
sx q[0];
rz(-0.17386359) q[0];
sx q[0];
rz(0.41282594) q[0];
rz(-1.6089571) q[1];
sx q[1];
rz(-2.2282579) q[1];
sx q[1];
rz(0.70721165) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4953046) q[0];
sx q[0];
rz(-2.2997059) q[0];
sx q[0];
rz(1.1524617) q[0];
x q[1];
rz(-0.74208547) q[2];
sx q[2];
rz(-1.4550503) q[2];
sx q[2];
rz(1.5884233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8726996) q[1];
sx q[1];
rz(-1.2352691) q[1];
sx q[1];
rz(-0.56443946) q[1];
rz(-pi) q[2];
rz(-1.3363289) q[3];
sx q[3];
rz(-2.5904561) q[3];
sx q[3];
rz(2.7296327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1327208) q[2];
sx q[2];
rz(-0.29942313) q[2];
sx q[2];
rz(-0.51187619) q[2];
rz(2.4608608) q[3];
sx q[3];
rz(-1.3635037) q[3];
sx q[3];
rz(-0.16630047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39636382) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(0.45528278) q[0];
rz(1.0633172) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(0.071970073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0531168) q[0];
sx q[0];
rz(-1.6040816) q[0];
sx q[0];
rz(-1.5831335) q[0];
rz(-pi) q[1];
x q[1];
rz(0.049111185) q[2];
sx q[2];
rz(-2.7002044) q[2];
sx q[2];
rz(2.9892444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1238511) q[1];
sx q[1];
rz(-1.2203274) q[1];
sx q[1];
rz(0.88092901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25983475) q[3];
sx q[3];
rz(-1.2502115) q[3];
sx q[3];
rz(-1.9966921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3528184) q[2];
sx q[2];
rz(-0.41663751) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(2.3173053) q[3];
sx q[3];
rz(-2.0485179) q[3];
sx q[3];
rz(0.85339671) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1085827) q[0];
sx q[0];
rz(-1.1065296) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(-1.4190326) q[1];
sx q[1];
rz(-0.483069) q[1];
sx q[1];
rz(-0.61666617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86980235) q[0];
sx q[0];
rz(-1.5951202) q[0];
sx q[0];
rz(-0.64016827) q[0];
rz(-1.9497112) q[2];
sx q[2];
rz(-0.99493631) q[2];
sx q[2];
rz(2.4873231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92041278) q[1];
sx q[1];
rz(-1.3592459) q[1];
sx q[1];
rz(-1.5736363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4343726) q[3];
sx q[3];
rz(-2.0248746) q[3];
sx q[3];
rz(-1.2421158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0509433) q[2];
sx q[2];
rz(-1.9885149) q[2];
sx q[2];
rz(2.0743745) q[2];
rz(2.0459335) q[3];
sx q[3];
rz(-1.2376384) q[3];
sx q[3];
rz(2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.511956) q[0];
sx q[0];
rz(-0.98573276) q[0];
sx q[0];
rz(-1.0712256) q[0];
rz(1.4818954) q[1];
sx q[1];
rz(-1.0922468) q[1];
sx q[1];
rz(-2.560871) q[1];
rz(-0.26255519) q[2];
sx q[2];
rz(-1.5468183) q[2];
sx q[2];
rz(0.47894947) q[2];
rz(1.7076013) q[3];
sx q[3];
rz(-1.8809263) q[3];
sx q[3];
rz(3.025007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
