OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8671829) q[0];
sx q[0];
rz(-2.6743968) q[0];
sx q[0];
rz(0.89611563) q[0];
rz(5.4652228) q[1];
sx q[1];
rz(4.735534) q[1];
sx q[1];
rz(8.4374333) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5817669) q[0];
sx q[0];
rz(-0.94115138) q[0];
sx q[0];
rz(2.6658076) q[0];
rz(-pi) q[1];
rz(2.5619123) q[2];
sx q[2];
rz(-0.90960056) q[2];
sx q[2];
rz(-1.5722317) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3442698) q[1];
sx q[1];
rz(-1.8472278) q[1];
sx q[1];
rz(1.396342) q[1];
rz(-1.7288858) q[3];
sx q[3];
rz(-1.2088547) q[3];
sx q[3];
rz(0.87495041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2200615) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(0.49911487) q[2];
rz(-2.3257997) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.65207425) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(-0.2336842) q[0];
rz(0.09985996) q[1];
sx q[1];
rz(-1.4870817) q[1];
sx q[1];
rz(-1.608009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.547116) q[0];
sx q[0];
rz(-2.9860382) q[0];
sx q[0];
rz(-2.2427325) q[0];
x q[1];
rz(-1.7560144) q[2];
sx q[2];
rz(-1.1389009) q[2];
sx q[2];
rz(0.76755953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74190226) q[1];
sx q[1];
rz(-0.61237017) q[1];
sx q[1];
rz(1.4111817) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30939524) q[3];
sx q[3];
rz(-1.1077048) q[3];
sx q[3];
rz(0.17406305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.097215501) q[2];
sx q[2];
rz(-2.2558687) q[2];
sx q[2];
rz(-0.054917939) q[2];
rz(0.6461668) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(0.072877876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570479) q[0];
sx q[0];
rz(-2.8962729) q[0];
sx q[0];
rz(-0.74485892) q[0];
rz(-1.2767731) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(0.863711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540032) q[0];
sx q[0];
rz(-1.1313725) q[0];
sx q[0];
rz(-1.5063398) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78869195) q[2];
sx q[2];
rz(-1.1967107) q[2];
sx q[2];
rz(-1.7497334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27093313) q[1];
sx q[1];
rz(-2.3133203) q[1];
sx q[1];
rz(0.67096295) q[1];
rz(-2.5287348) q[3];
sx q[3];
rz(-0.75988673) q[3];
sx q[3];
rz(-2.5155544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6700217) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(-0.395533) q[2];
rz(1.0223355) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763497) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(-3.0856207) q[0];
rz(2.5296027) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(-0.8459808) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6406785) q[0];
sx q[0];
rz(-2.0862824) q[0];
sx q[0];
rz(2.8436996) q[0];
x q[1];
rz(0.34807713) q[2];
sx q[2];
rz(-0.98053369) q[2];
sx q[2];
rz(-1.2952309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27056672) q[1];
sx q[1];
rz(-1.881384) q[1];
sx q[1];
rz(1.3695716) q[1];
rz(-pi) q[2];
rz(1.0958798) q[3];
sx q[3];
rz(-1.2327415) q[3];
sx q[3];
rz(-2.2477704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(-2.4212627) q[2];
rz(2.7885041) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(-2.8928355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40097749) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(2.154696) q[0];
rz(1.1445716) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(-2.1867337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3026149) q[0];
sx q[0];
rz(-2.2377126) q[0];
sx q[0];
rz(3.1083005) q[0];
rz(-pi) q[1];
rz(-1.0748765) q[2];
sx q[2];
rz(-0.62037797) q[2];
sx q[2];
rz(1.4324783) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2209619) q[1];
sx q[1];
rz(-0.60505962) q[1];
sx q[1];
rz(-0.69843881) q[1];
rz(-pi) q[2];
rz(-0.14820672) q[3];
sx q[3];
rz(-1.3199521) q[3];
sx q[3];
rz(0.36484066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6285051) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(-1.0465735) q[2];
rz(-2.4322677) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76383048) q[0];
sx q[0];
rz(-1.4077633) q[0];
sx q[0];
rz(-2.8705226) q[0];
rz(0.079004869) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(-1.7105506) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56667152) q[0];
sx q[0];
rz(-1.780176) q[0];
sx q[0];
rz(1.4230595) q[0];
rz(-pi) q[1];
rz(2.3777804) q[2];
sx q[2];
rz(-2.4226563) q[2];
sx q[2];
rz(0.97455183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5137427) q[1];
sx q[1];
rz(-1.8166891) q[1];
sx q[1];
rz(-2.4469083) q[1];
rz(1.920041) q[3];
sx q[3];
rz(-2.6995097) q[3];
sx q[3];
rz(0.52607049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0926823) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(3.1177706) q[2];
rz(1.228099) q[3];
sx q[3];
rz(-0.3796328) q[3];
sx q[3];
rz(-2.9983799) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.989885) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(-0.09224961) q[0];
rz(-0.30091885) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-1.0801962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5841519) q[0];
sx q[0];
rz(-1.6655465) q[0];
sx q[0];
rz(1.9173724) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.779498) q[2];
sx q[2];
rz(-1.9685192) q[2];
sx q[2];
rz(2.5366572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7799973) q[1];
sx q[1];
rz(-1.4228586) q[1];
sx q[1];
rz(1.5725488) q[1];
rz(0.37410847) q[3];
sx q[3];
rz(-2.4112933) q[3];
sx q[3];
rz(2.7299943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-2.7907382) q[2];
sx q[2];
rz(0.62823137) q[2];
rz(-2.175711) q[3];
sx q[3];
rz(-1.5768257) q[3];
sx q[3];
rz(-2.8204744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0308762) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(-1.734717) q[0];
rz(0.12786099) q[1];
sx q[1];
rz(-1.4297994) q[1];
sx q[1];
rz(2.0157287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073407452) q[0];
sx q[0];
rz(-0.74755423) q[0];
sx q[0];
rz(0.27184591) q[0];
rz(-pi) q[1];
rz(-2.0225384) q[2];
sx q[2];
rz(-1.5875419) q[2];
sx q[2];
rz(-2.708205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5347157) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(1.5061395) q[1];
rz(-pi) q[2];
rz(0.58784318) q[3];
sx q[3];
rz(-2.4333262) q[3];
sx q[3];
rz(1.4204587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0109978) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(-3.0478743) q[2];
rz(1.7081918) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4505287) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(-1.0040671) q[0];
rz(2.0380691) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(-1.4080661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44847261) q[0];
sx q[0];
rz(-1.5647443) q[0];
sx q[0];
rz(0.5151185) q[0];
x q[1];
rz(0.33317487) q[2];
sx q[2];
rz(-2.8123724) q[2];
sx q[2];
rz(-1.7235989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50251657) q[1];
sx q[1];
rz(-2.6682819) q[1];
sx q[1];
rz(1.972354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.086661913) q[3];
sx q[3];
rz(-2.6209957) q[3];
sx q[3];
rz(-2.749032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7496926) q[2];
sx q[2];
rz(-1.6204648) q[2];
sx q[2];
rz(-2.5938972) q[2];
rz(2.8294166) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(-1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80426973) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(-0.7789337) q[0];
rz(0.3745105) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(-2.3096854) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9269126) q[0];
sx q[0];
rz(-2.2601366) q[0];
sx q[0];
rz(-1.229152) q[0];
rz(-pi) q[1];
rz(0.7829297) q[2];
sx q[2];
rz(-2.1179869) q[2];
sx q[2];
rz(1.3732571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0613855) q[1];
sx q[1];
rz(-2.204772) q[1];
sx q[1];
rz(-0.69532077) q[1];
rz(-pi) q[2];
rz(-0.030422525) q[3];
sx q[3];
rz(-0.796954) q[3];
sx q[3];
rz(1.6773633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0348009) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(-0.0084776004) q[2];
rz(2.7360385) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(-1.6976374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99210284) q[0];
sx q[0];
rz(-1.6076417) q[0];
sx q[0];
rz(-2.3544307) q[0];
rz(-1.0943195) q[1];
sx q[1];
rz(-0.37847395) q[1];
sx q[1];
rz(2.7056221) q[1];
rz(-3.0158184) q[2];
sx q[2];
rz(-1.3604506) q[2];
sx q[2];
rz(2.4764555) q[2];
rz(1.5658436) q[3];
sx q[3];
rz(-1.4765783) q[3];
sx q[3];
rz(0.12259132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
