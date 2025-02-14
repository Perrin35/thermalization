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
rz(0.12198099) q[0];
sx q[0];
rz(-0.044959083) q[0];
sx q[0];
rz(-3.0461351) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(3.0756693) q[1];
sx q[1];
rz(8.8311721) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55768591) q[0];
sx q[0];
rz(-0.61827165) q[0];
sx q[0];
rz(2.4960211) q[0];
rz(-pi) q[1];
rz(0.66887899) q[2];
sx q[2];
rz(-2.9496585) q[2];
sx q[2];
rz(-0.30489433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7194743) q[1];
sx q[1];
rz(-2.1917043) q[1];
sx q[1];
rz(1.3139753) q[1];
rz(2.8659561) q[3];
sx q[3];
rz(-1.0616117) q[3];
sx q[3];
rz(-1.0643626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.060071271) q[2];
sx q[2];
rz(-1.2520484) q[2];
sx q[2];
rz(0.56212765) q[2];
rz(2.8020322) q[3];
sx q[3];
rz(-1.6189251) q[3];
sx q[3];
rz(-1.08574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.958441) q[0];
sx q[0];
rz(-0.14587942) q[0];
sx q[0];
rz(-1.1613783) q[0];
rz(2.9005652) q[1];
sx q[1];
rz(-2.3190505) q[1];
sx q[1];
rz(2.8384812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35251543) q[0];
sx q[0];
rz(-0.20931986) q[0];
sx q[0];
rz(1.2417488) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70485665) q[2];
sx q[2];
rz(-2.9508756) q[2];
sx q[2];
rz(-2.5097367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1352586) q[1];
sx q[1];
rz(-2.5016682) q[1];
sx q[1];
rz(2.695822) q[1];
rz(-pi) q[2];
rz(2.3645401) q[3];
sx q[3];
rz(-2.5133479) q[3];
sx q[3];
rz(-0.45690445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8384398) q[2];
sx q[2];
rz(-1.2829605) q[2];
sx q[2];
rz(-2.8008833) q[2];
rz(-0.43863145) q[3];
sx q[3];
rz(-2.3352968) q[3];
sx q[3];
rz(3.0830234) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71857798) q[0];
sx q[0];
rz(-2.5231762) q[0];
sx q[0];
rz(-2.3037236) q[0];
rz(-0.92309976) q[1];
sx q[1];
rz(-1.2078614) q[1];
sx q[1];
rz(2.7941678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14590696) q[0];
sx q[0];
rz(-1.4769698) q[0];
sx q[0];
rz(-0.037414649) q[0];
rz(1.8272039) q[2];
sx q[2];
rz(-1.5359582) q[2];
sx q[2];
rz(0.23171356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4975884) q[1];
sx q[1];
rz(-2.1198744) q[1];
sx q[1];
rz(1.0058897) q[1];
rz(-pi) q[2];
rz(2.8108571) q[3];
sx q[3];
rz(-1.8285255) q[3];
sx q[3];
rz(2.7455519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0190987) q[2];
sx q[2];
rz(-3.0235897) q[2];
sx q[2];
rz(-0.72354358) q[2];
rz(-1.8540234) q[3];
sx q[3];
rz(-0.55793327) q[3];
sx q[3];
rz(3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36640722) q[0];
sx q[0];
rz(-2.7025096) q[0];
sx q[0];
rz(-2.0570237) q[0];
rz(1.9547801) q[1];
sx q[1];
rz(-1.908952) q[1];
sx q[1];
rz(1.67217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0908244) q[0];
sx q[0];
rz(-1.6009129) q[0];
sx q[0];
rz(-1.5876905) q[0];
rz(-pi) q[1];
rz(-2.3575063) q[2];
sx q[2];
rz(-0.52496451) q[2];
sx q[2];
rz(2.3388179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.332843) q[1];
sx q[1];
rz(-1.9522138) q[1];
sx q[1];
rz(-1.0709958) q[1];
x q[2];
rz(0.74923781) q[3];
sx q[3];
rz(-2.6686274) q[3];
sx q[3];
rz(1.7974896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12545086) q[2];
sx q[2];
rz(-0.90105385) q[2];
sx q[2];
rz(-0.93552843) q[2];
rz(2.9901796) q[3];
sx q[3];
rz(-2.1617523) q[3];
sx q[3];
rz(2.6305731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147515) q[0];
sx q[0];
rz(-1.5578569) q[0];
sx q[0];
rz(0.75411183) q[0];
rz(0.61112815) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(0.67620826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5181715) q[0];
sx q[0];
rz(-1.8966798) q[0];
sx q[0];
rz(-0.17101479) q[0];
x q[1];
rz(1.0173747) q[2];
sx q[2];
rz(-0.62741919) q[2];
sx q[2];
rz(1.8946033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0724746) q[1];
sx q[1];
rz(-2.9058911) q[1];
sx q[1];
rz(0.75171296) q[1];
rz(-pi) q[2];
rz(1.6531472) q[3];
sx q[3];
rz(-2.1550278) q[3];
sx q[3];
rz(2.1305934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0040697441) q[2];
sx q[2];
rz(-2.4293032) q[2];
sx q[2];
rz(2.3884657) q[2];
rz(-2.5774041) q[3];
sx q[3];
rz(-2.0337532) q[3];
sx q[3];
rz(-2.3851725) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18993987) q[0];
sx q[0];
rz(-0.48007444) q[0];
sx q[0];
rz(-0.22115627) q[0];
rz(3.0033374) q[1];
sx q[1];
rz(-1.4555376) q[1];
sx q[1];
rz(-2.5855605) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0079239844) q[0];
sx q[0];
rz(-1.3268262) q[0];
sx q[0];
rz(-1.7708112) q[0];
rz(-pi) q[1];
rz(-1.1821317) q[2];
sx q[2];
rz(-1.2209397) q[2];
sx q[2];
rz(-1.4887202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41084633) q[1];
sx q[1];
rz(-1.7241983) q[1];
sx q[1];
rz(-1.549841) q[1];
x q[2];
rz(2.0222072) q[3];
sx q[3];
rz(-1.4638293) q[3];
sx q[3];
rz(0.2052923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14744645) q[2];
sx q[2];
rz(-0.7615971) q[2];
sx q[2];
rz(-2.2046294) q[2];
rz(2.7312036) q[3];
sx q[3];
rz(-2.161721) q[3];
sx q[3];
rz(-2.4858937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74528247) q[0];
sx q[0];
rz(-0.87438011) q[0];
sx q[0];
rz(2.1790047) q[0];
rz(2.4707322) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(0.76505351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214692) q[0];
sx q[0];
rz(-1.7304106) q[0];
sx q[0];
rz(-1.1667211) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1704438) q[2];
sx q[2];
rz(-0.12242386) q[2];
sx q[2];
rz(-2.3442307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1598164) q[1];
sx q[1];
rz(-0.66513956) q[1];
sx q[1];
rz(-0.90250166) q[1];
x q[2];
rz(1.2575862) q[3];
sx q[3];
rz(-0.38496415) q[3];
sx q[3];
rz(2.0244618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44136167) q[2];
sx q[2];
rz(-2.2141778) q[2];
sx q[2];
rz(0.78511635) q[2];
rz(-2.650812) q[3];
sx q[3];
rz(-1.9793341) q[3];
sx q[3];
rz(2.6665915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15071507) q[0];
sx q[0];
rz(-0.10902037) q[0];
sx q[0];
rz(0.14608598) q[0];
rz(-0.27169216) q[1];
sx q[1];
rz(-0.79333317) q[1];
sx q[1];
rz(-2.9300516) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95446996) q[0];
sx q[0];
rz(-0.7899219) q[0];
sx q[0];
rz(0.23440897) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8604989) q[2];
sx q[2];
rz(-0.38662505) q[2];
sx q[2];
rz(-2.7747216) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.817343) q[1];
sx q[1];
rz(-1.5012263) q[1];
sx q[1];
rz(-2.9051498) q[1];
rz(-pi) q[2];
rz(1.7320427) q[3];
sx q[3];
rz(-2.5354156) q[3];
sx q[3];
rz(2.2796749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81296706) q[2];
sx q[2];
rz(-0.93572891) q[2];
sx q[2];
rz(-2.3198371) q[2];
rz(1.8030608) q[3];
sx q[3];
rz(-1.0889564) q[3];
sx q[3];
rz(0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5271135) q[0];
sx q[0];
rz(-0.69632691) q[0];
sx q[0];
rz(-1.7823727) q[0];
rz(0.95787734) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(2.8639796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045317) q[0];
sx q[0];
rz(-0.87420428) q[0];
sx q[0];
rz(0.84203984) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3657655) q[2];
sx q[2];
rz(-1.7971353) q[2];
sx q[2];
rz(-0.090858484) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0467068) q[1];
sx q[1];
rz(-1.1562655) q[1];
sx q[1];
rz(1.4317896) q[1];
rz(-1.3986392) q[3];
sx q[3];
rz(-1.1973698) q[3];
sx q[3];
rz(0.40478727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9684888) q[2];
sx q[2];
rz(-0.59440458) q[2];
sx q[2];
rz(-1.8572726) q[2];
rz(0.77749085) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(2.8418181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023329968) q[0];
sx q[0];
rz(-1.6602004) q[0];
sx q[0];
rz(3.135664) q[0];
rz(-1.802035) q[1];
sx q[1];
rz(-1.0362933) q[1];
sx q[1];
rz(2.6609227) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0927785) q[0];
sx q[0];
rz(-1.5869821) q[0];
sx q[0];
rz(-1.556514) q[0];
rz(-pi) q[1];
rz(-2.2193682) q[2];
sx q[2];
rz(-1.2184452) q[2];
sx q[2];
rz(1.00204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.103616) q[1];
sx q[1];
rz(-2.3876752) q[1];
sx q[1];
rz(-0.62289728) q[1];
x q[2];
rz(-2.0391614) q[3];
sx q[3];
rz(-2.8448555) q[3];
sx q[3];
rz(0.67978978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9504451) q[2];
sx q[2];
rz(-0.55926776) q[2];
sx q[2];
rz(-0.17678235) q[2];
rz(-2.2117129) q[3];
sx q[3];
rz(-1.1888489) q[3];
sx q[3];
rz(2.1117579) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4491691) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(-2.069166) q[1];
sx q[1];
rz(-1.4785531) q[1];
sx q[1];
rz(-1.4428152) q[1];
rz(1.564413) q[2];
sx q[2];
rz(-1.0789568) q[2];
sx q[2];
rz(3.073624) q[2];
rz(-1.4789875) q[3];
sx q[3];
rz(-0.96755618) q[3];
sx q[3];
rz(0.70152828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
