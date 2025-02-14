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
rz(1.0072768) q[0];
sx q[0];
rz(-0.64581031) q[0];
sx q[0];
rz(0.3663775) q[0];
rz(0.74908787) q[1];
sx q[1];
rz(4.7685342) q[1];
sx q[1];
rz(9.2821477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86038268) q[0];
sx q[0];
rz(-1.027635) q[0];
sx q[0];
rz(-2.5090911) q[0];
x q[1];
rz(-2.9171474) q[2];
sx q[2];
rz(-1.7068752) q[2];
sx q[2];
rz(-0.70101695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4971338) q[1];
sx q[1];
rz(-2.6503149) q[1];
sx q[1];
rz(1.8257479) q[1];
rz(1.8964389) q[3];
sx q[3];
rz(-2.0228141) q[3];
sx q[3];
rz(-2.8309517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7819405) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(-1.5521607) q[2];
rz(-1.747067) q[3];
sx q[3];
rz(-2.7701869) q[3];
sx q[3];
rz(2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29779103) q[0];
sx q[0];
rz(-2.5318635) q[0];
sx q[0];
rz(-2.8363256) q[0];
rz(1.3097395) q[1];
sx q[1];
rz(-0.42116183) q[1];
sx q[1];
rz(3.0895244) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48114313) q[0];
sx q[0];
rz(-1.423712) q[0];
sx q[0];
rz(1.4597465) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7676653) q[2];
sx q[2];
rz(-1.2902593) q[2];
sx q[2];
rz(-0.032023059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4365719) q[1];
sx q[1];
rz(-1.8364882) q[1];
sx q[1];
rz(0.17036856) q[1];
rz(-pi) q[2];
rz(1.4036739) q[3];
sx q[3];
rz(-0.67490904) q[3];
sx q[3];
rz(3.1233146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0569968) q[2];
sx q[2];
rz(-1.8085542) q[2];
sx q[2];
rz(1.7219273) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(1.4066633) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97515714) q[0];
sx q[0];
rz(-1.1844013) q[0];
sx q[0];
rz(1.3982406) q[0];
rz(0.94683975) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(1.2907226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7311659) q[0];
sx q[0];
rz(-1.7058347) q[0];
sx q[0];
rz(0.076083994) q[0];
x q[1];
rz(2.0507347) q[2];
sx q[2];
rz(-2.2971114) q[2];
sx q[2];
rz(-1.2588071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5980893) q[1];
sx q[1];
rz(-3.0687947) q[1];
sx q[1];
rz(-1.9803712) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3685143) q[3];
sx q[3];
rz(-2.0401388) q[3];
sx q[3];
rz(-2.3465867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1281841) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(0.95995963) q[2];
rz(0.32402447) q[3];
sx q[3];
rz(-0.5158546) q[3];
sx q[3];
rz(0.44017756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1337166) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(-2.9414951) q[0];
rz(0.61087418) q[1];
sx q[1];
rz(-2.4433544) q[1];
sx q[1];
rz(-2.3470338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66629529) q[0];
sx q[0];
rz(-1.4213511) q[0];
sx q[0];
rz(1.4578865) q[0];
rz(2.9616293) q[2];
sx q[2];
rz(-0.90562253) q[2];
sx q[2];
rz(0.53097224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9640089) q[1];
sx q[1];
rz(-1.4191886) q[1];
sx q[1];
rz(1.2535918) q[1];
x q[2];
rz(0.67312397) q[3];
sx q[3];
rz(-0.84852058) q[3];
sx q[3];
rz(-0.22658843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8841298) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(1.3563159) q[3];
sx q[3];
rz(-1.9808199) q[3];
sx q[3];
rz(0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2785579) q[0];
sx q[0];
rz(-2.886241) q[0];
sx q[0];
rz(-0.72364664) q[0];
rz(-0.10069314) q[1];
sx q[1];
rz(-2.0932525) q[1];
sx q[1];
rz(0.06631276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6479939) q[0];
sx q[0];
rz(-1.0033404) q[0];
sx q[0];
rz(-3.0667414) q[0];
rz(0.91513779) q[2];
sx q[2];
rz(-1.1205289) q[2];
sx q[2];
rz(1.428626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1182013) q[1];
sx q[1];
rz(-1.4033591) q[1];
sx q[1];
rz(-0.91806478) q[1];
x q[2];
rz(1.2454493) q[3];
sx q[3];
rz(-1.950437) q[3];
sx q[3];
rz(1.0285447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.940332) q[2];
sx q[2];
rz(-1.9278434) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(-0.39441937) q[3];
sx q[3];
rz(-0.68587488) q[3];
sx q[3];
rz(-2.7421537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0938996) q[0];
sx q[0];
rz(-2.221929) q[0];
sx q[0];
rz(-3.1374875) q[0];
rz(-2.507569) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(0.9563458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8656971) q[0];
sx q[0];
rz(-1.529954) q[0];
sx q[0];
rz(2.9581822) q[0];
x q[1];
rz(1.1460828) q[2];
sx q[2];
rz(-2.0419952) q[2];
sx q[2];
rz(1.4482519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8022393) q[1];
sx q[1];
rz(-1.5763073) q[1];
sx q[1];
rz(0.27779419) q[1];
rz(-pi) q[2];
rz(0.31305571) q[3];
sx q[3];
rz(-2.7767973) q[3];
sx q[3];
rz(-0.49653253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6213106) q[2];
sx q[2];
rz(-1.7063528) q[2];
sx q[2];
rz(-0.83462805) q[2];
rz(2.7708715) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(0.71681517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1762367) q[0];
sx q[0];
rz(-3.1270202) q[0];
sx q[0];
rz(0.45005774) q[0];
rz(0.4484446) q[1];
sx q[1];
rz(-2.3095755) q[1];
sx q[1];
rz(-2.4041596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56313228) q[0];
sx q[0];
rz(-1.9942584) q[0];
sx q[0];
rz(2.9388301) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.040410553) q[2];
sx q[2];
rz(-1.4295661) q[2];
sx q[2];
rz(-0.34858957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1597643) q[1];
sx q[1];
rz(-0.48063403) q[1];
sx q[1];
rz(-2.3932061) q[1];
rz(-1.3025018) q[3];
sx q[3];
rz(-2.3327069) q[3];
sx q[3];
rz(-0.20352645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8488778) q[2];
sx q[2];
rz(-0.60056168) q[2];
sx q[2];
rz(-0.71504492) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.4628937) q[3];
sx q[3];
rz(-1.3983294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364863) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(-0.35838321) q[0];
rz(0.37250039) q[1];
sx q[1];
rz(-1.3595711) q[1];
sx q[1];
rz(0.43982664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7067703) q[0];
sx q[0];
rz(-0.79688886) q[0];
sx q[0];
rz(-3.0838941) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.110564) q[2];
sx q[2];
rz(-1.2046736) q[2];
sx q[2];
rz(-1.5535959) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.81573) q[1];
sx q[1];
rz(-1.7472032) q[1];
sx q[1];
rz(0.26858624) q[1];
rz(-2.8183455) q[3];
sx q[3];
rz(-1.6222937) q[3];
sx q[3];
rz(-0.53158376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0899352) q[2];
sx q[2];
rz(-0.097044162) q[2];
sx q[2];
rz(0.70866054) q[2];
rz(2.8900201) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(-3.1075509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8468903) q[0];
sx q[0];
rz(-2.1069694) q[0];
sx q[0];
rz(2.5501472) q[0];
rz(-0.017631831) q[1];
sx q[1];
rz(-1.0523187) q[1];
sx q[1];
rz(-0.10094053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4288947) q[0];
sx q[0];
rz(-2.7980045) q[0];
sx q[0];
rz(1.1593007) q[0];
rz(1.5315311) q[2];
sx q[2];
rz(-2.2071243) q[2];
sx q[2];
rz(-2.1992342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7314607) q[1];
sx q[1];
rz(-1.6900151) q[1];
sx q[1];
rz(1.340926) q[1];
rz(-pi) q[2];
rz(0.89102192) q[3];
sx q[3];
rz(-1.4263065) q[3];
sx q[3];
rz(1.4308892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4883604) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(1.4867268) q[2];
rz(2.9871353) q[3];
sx q[3];
rz(-2.0363225) q[3];
sx q[3];
rz(0.36623335) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79009295) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(-1.0737786) q[0];
rz(2.0692661) q[1];
sx q[1];
rz(-0.25389478) q[1];
sx q[1];
rz(-0.059018746) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71145161) q[0];
sx q[0];
rz(-1.7699598) q[0];
sx q[0];
rz(1.6047305) q[0];
rz(-pi) q[1];
rz(-2.8779936) q[2];
sx q[2];
rz(-1.0754183) q[2];
sx q[2];
rz(-1.6558629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37009753) q[1];
sx q[1];
rz(-2.1019249) q[1];
sx q[1];
rz(-0.70651502) q[1];
rz(-pi) q[2];
x q[2];
rz(2.288184) q[3];
sx q[3];
rz(-2.0730258) q[3];
sx q[3];
rz(-2.2706007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3286288) q[2];
sx q[2];
rz(-1.8495411) q[2];
sx q[2];
rz(0.22895075) q[2];
rz(1.1936584) q[3];
sx q[3];
rz(-0.3242068) q[3];
sx q[3];
rz(1.4540023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1345632) q[0];
sx q[0];
rz(-0.80637359) q[0];
sx q[0];
rz(-0.73643186) q[0];
rz(1.7560584) q[1];
sx q[1];
rz(-1.3722739) q[1];
sx q[1];
rz(-1.3442232) q[1];
rz(-0.028111721) q[2];
sx q[2];
rz(-2.1098469) q[2];
sx q[2];
rz(-2.5218365) q[2];
rz(-1.368461) q[3];
sx q[3];
rz(-2.1740395) q[3];
sx q[3];
rz(-2.5185829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
