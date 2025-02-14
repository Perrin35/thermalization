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
rz(0.8402549) q[0];
sx q[0];
rz(4.1794887) q[0];
sx q[0];
rz(10.269796) q[0];
rz(-2.4294699) q[1];
sx q[1];
rz(-1.0016088) q[1];
sx q[1];
rz(1.6460302) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9121542) q[0];
sx q[0];
rz(-2.144838) q[0];
sx q[0];
rz(0.50749929) q[0];
rz(-pi) q[1];
rz(1.4207815) q[2];
sx q[2];
rz(-1.1910465) q[2];
sx q[2];
rz(-0.94204547) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.597376) q[1];
sx q[1];
rz(-2.4565182) q[1];
sx q[1];
rz(2.3822576) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2949307) q[3];
sx q[3];
rz(-1.1893236) q[3];
sx q[3];
rz(-1.4700397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40319765) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(0.5564059) q[2];
rz(0.49499908) q[3];
sx q[3];
rz(-0.35111108) q[3];
sx q[3];
rz(-0.88465148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2719088) q[0];
sx q[0];
rz(-3.0392635) q[0];
sx q[0];
rz(2.5129357) q[0];
rz(2.2643845) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(2.8884851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0639187) q[0];
sx q[0];
rz(-1.5787933) q[0];
sx q[0];
rz(-1.5551989) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2727401) q[2];
sx q[2];
rz(-1.9949081) q[2];
sx q[2];
rz(-1.9587245) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4483466) q[1];
sx q[1];
rz(-2.9102059) q[1];
sx q[1];
rz(0.53129249) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3023002) q[3];
sx q[3];
rz(-1.3802234) q[3];
sx q[3];
rz(0.96123248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6796598) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(1.0665464) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(-2.2333142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.2147373) q[0];
sx q[0];
rz(-0.94811386) q[0];
sx q[0];
rz(-2.1562449) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-0.99291283) q[1];
sx q[1];
rz(-2.6148112) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5827127) q[0];
sx q[0];
rz(-2.0979019) q[0];
sx q[0];
rz(-0.096313535) q[0];
rz(2.5622732) q[2];
sx q[2];
rz(-1.123482) q[2];
sx q[2];
rz(1.2944702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0635447) q[1];
sx q[1];
rz(-1.9920262) q[1];
sx q[1];
rz(2.9560412) q[1];
rz(-0.95895264) q[3];
sx q[3];
rz(-2.4482895) q[3];
sx q[3];
rz(2.7075651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3452611) q[2];
sx q[2];
rz(-2.3806206) q[2];
sx q[2];
rz(2.995028) q[2];
rz(0.86756271) q[3];
sx q[3];
rz(-0.87351322) q[3];
sx q[3];
rz(-0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.72916156) q[0];
sx q[0];
rz(-2.6456092) q[0];
sx q[0];
rz(-0.55369401) q[0];
rz(1.4750534) q[1];
sx q[1];
rz(-2.8429884) q[1];
sx q[1];
rz(-0.24436229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840773) q[0];
sx q[0];
rz(-1.2173442) q[0];
sx q[0];
rz(1.9259324) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9072794) q[2];
sx q[2];
rz(-2.3762616) q[2];
sx q[2];
rz(-2.825277) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3267388) q[1];
sx q[1];
rz(-2.5757901) q[1];
sx q[1];
rz(-0.71235719) q[1];
rz(-3.1126106) q[3];
sx q[3];
rz(-0.71235114) q[3];
sx q[3];
rz(2.1952352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.060140572) q[2];
sx q[2];
rz(-2.0144561) q[2];
sx q[2];
rz(-0.78090182) q[2];
rz(2.7447356) q[3];
sx q[3];
rz(-0.01440993) q[3];
sx q[3];
rz(-2.0675596) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395372) q[0];
sx q[0];
rz(-0.1665512) q[0];
sx q[0];
rz(-0.2567513) q[0];
rz(0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(0.94388747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092800822) q[0];
sx q[0];
rz(-1.208713) q[0];
sx q[0];
rz(1.0975773) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19380865) q[2];
sx q[2];
rz(-1.8558981) q[2];
sx q[2];
rz(1.5616035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5748716) q[1];
sx q[1];
rz(-1.566426) q[1];
sx q[1];
rz(2.6019051) q[1];
rz(2.6750608) q[3];
sx q[3];
rz(-2.0781019) q[3];
sx q[3];
rz(-3.0477085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93861598) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(-0.20982783) q[2];
rz(-0.046791568) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-2.7941424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.2030075) q[0];
rz(1.5164392) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(0.032940544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0896785) q[0];
sx q[0];
rz(-2.5915049) q[0];
sx q[0];
rz(2.1142531) q[0];
x q[1];
rz(-2.2896272) q[2];
sx q[2];
rz(-2.4174066) q[2];
sx q[2];
rz(0.23978309) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60073416) q[1];
sx q[1];
rz(-1.2493397) q[1];
sx q[1];
rz(-0.58166166) q[1];
rz(0.83077151) q[3];
sx q[3];
rz(-2.3099595) q[3];
sx q[3];
rz(-2.3297854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5861627) q[2];
sx q[2];
rz(-2.1383492) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(0.14719506) q[3];
sx q[3];
rz(-2.9087524) q[3];
sx q[3];
rz(1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066605695) q[0];
sx q[0];
rz(-2.7050278) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(-2.1287411) q[1];
sx q[1];
rz(-1.9722936) q[1];
sx q[1];
rz(0.036651932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3507639) q[0];
sx q[0];
rz(-1.1393271) q[0];
sx q[0];
rz(0.58833265) q[0];
rz(1.4674856) q[2];
sx q[2];
rz(-1.9983618) q[2];
sx q[2];
rz(0.72139886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0554296) q[1];
sx q[1];
rz(-0.87172548) q[1];
sx q[1];
rz(0.12075874) q[1];
rz(0.030825214) q[3];
sx q[3];
rz(-1.5753058) q[3];
sx q[3];
rz(2.4140178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6645633) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-0.089740962) q[2];
rz(1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(1.4242273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.4795714) q[0];
sx q[0];
rz(-2.5630072) q[0];
sx q[0];
rz(1.660996) q[0];
rz(-1.4501694) q[1];
sx q[1];
rz(-0.65933508) q[1];
sx q[1];
rz(2.3816542) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1131033) q[0];
sx q[0];
rz(-0.11664243) q[0];
sx q[0];
rz(-1.4820497) q[0];
rz(3.0133661) q[2];
sx q[2];
rz(-2.5094446) q[2];
sx q[2];
rz(-1.3299483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41553283) q[1];
sx q[1];
rz(-0.59846717) q[1];
sx q[1];
rz(2.3319671) q[1];
x q[2];
rz(1.8283056) q[3];
sx q[3];
rz(-1.0777359) q[3];
sx q[3];
rz(0.3230394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8370168) q[2];
sx q[2];
rz(-1.2695856) q[2];
sx q[2];
rz(-0.059583511) q[2];
rz(2.7130821) q[3];
sx q[3];
rz(-2.2969552) q[3];
sx q[3];
rz(-0.60215157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483599) q[0];
sx q[0];
rz(-2.8587274) q[0];
sx q[0];
rz(-0.34348139) q[0];
rz(2.9454625) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(0.79998618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50757396) q[0];
sx q[0];
rz(-0.91092604) q[0];
sx q[0];
rz(2.6150661) q[0];
rz(-0.38651379) q[2];
sx q[2];
rz(-0.59745759) q[2];
sx q[2];
rz(0.51624417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5243036) q[1];
sx q[1];
rz(-2.5802748) q[1];
sx q[1];
rz(0.08331068) q[1];
rz(-3.0723582) q[3];
sx q[3];
rz(-2.2524009) q[3];
sx q[3];
rz(2.414976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34667748) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(-0.22870341) q[2];
rz(1.4165357) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(-0.67030877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59920853) q[0];
sx q[0];
rz(-0.51775652) q[0];
sx q[0];
rz(1.1826578) q[0];
rz(0.54284894) q[1];
sx q[1];
rz(-0.12195568) q[1];
sx q[1];
rz(0.45546946) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96129721) q[0];
sx q[0];
rz(-2.0019128) q[0];
sx q[0];
rz(-0.53133734) q[0];
rz(1.1588016) q[2];
sx q[2];
rz(-2.6557026) q[2];
sx q[2];
rz(-2.0249572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7389981) q[1];
sx q[1];
rz(-1.4345503) q[1];
sx q[1];
rz(-0.46936492) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72584589) q[3];
sx q[3];
rz(-2.5739087) q[3];
sx q[3];
rz(-0.7782225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9246284) q[2];
sx q[2];
rz(-0.89322007) q[2];
sx q[2];
rz(-2.7157937) q[2];
rz(-0.18481542) q[3];
sx q[3];
rz(-2.5033689) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.4375147) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(-2.6620445) q[1];
sx q[1];
rz(-1.5794812) q[1];
sx q[1];
rz(-1.5522122) q[1];
rz(-1.7935971) q[2];
sx q[2];
rz(-2.1542414) q[2];
sx q[2];
rz(0.20050394) q[2];
rz(-0.14397022) q[3];
sx q[3];
rz(-2.1767749) q[3];
sx q[3];
rz(-0.53713606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
