OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0159338) q[0];
sx q[0];
rz(-2.4056245) q[0];
sx q[0];
rz(-2.486881) q[0];
rz(2.4891698) q[2];
sx q[2];
rz(-2.4180275) q[2];
sx q[2];
rz(-0.10318081) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9846676) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.9763293) q[1];
rz(-1.3017544) q[3];
sx q[3];
rz(-1.6645414) q[3];
sx q[3];
rz(0.71414381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(0.75749767) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(-2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(-0.27145162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4685681) q[0];
sx q[0];
rz(-2.4789171) q[0];
sx q[0];
rz(-2.2492692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2618622) q[2];
sx q[2];
rz(-2.0663107) q[2];
sx q[2];
rz(-2.3625284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8703692) q[1];
sx q[1];
rz(-1.1961812) q[1];
sx q[1];
rz(-0.80925525) q[1];
rz(1.9430964) q[3];
sx q[3];
rz(-2.3447403) q[3];
sx q[3];
rz(-0.70355319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(2.2180637) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.7318168) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554879) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(-0.79746042) q[0];
rz(-pi) q[1];
rz(0.44720165) q[2];
sx q[2];
rz(-0.7913835) q[2];
sx q[2];
rz(2.7176822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8612954) q[1];
sx q[1];
rz(-2.0164844) q[1];
sx q[1];
rz(-0.46511005) q[1];
rz(0.61061065) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(-3.1162457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(-0.96281111) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(3.1406291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9199333) q[0];
sx q[0];
rz(-0.98826212) q[0];
sx q[0];
rz(2.4525053) q[0];
rz(-pi) q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5713072) q[2];
sx q[2];
rz(2.309547) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2767267) q[1];
sx q[1];
rz(-2.556567) q[1];
sx q[1];
rz(-1.8086955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6603052) q[3];
sx q[3];
rz(-2.138391) q[3];
sx q[3];
rz(-1.2168509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-0.11165079) q[2];
rz(-2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(0.26062632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5379932) q[0];
sx q[0];
rz(-1.4289745) q[0];
sx q[0];
rz(1.5285368) q[0];
x q[1];
rz(-1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(-1.7340811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6071636) q[1];
sx q[1];
rz(-0.57463127) q[1];
sx q[1];
rz(-0.032392153) q[1];
rz(-pi) q[2];
rz(0.77633206) q[3];
sx q[3];
rz(-2.9388802) q[3];
sx q[3];
rz(2.4622038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2300718) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(-1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.7274436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0522239) q[0];
sx q[0];
rz(-1.5234689) q[0];
sx q[0];
rz(-2.8663859) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44712375) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(-2.2396357) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8687739) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(0.52656071) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9023444) q[3];
sx q[3];
rz(-2.9132531) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(-3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(1.3940575) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(2.7391403) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50169045) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-2.7946266) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(2.7260821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35343364) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-0.56281705) q[1];
rz(-pi) q[2];
rz(-1.2225371) q[3];
sx q[3];
rz(-1.5969443) q[3];
sx q[3];
rz(0.98046434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1272614) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-1.0345116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1326133) q[0];
sx q[0];
rz(-1.4011369) q[0];
sx q[0];
rz(-1.8011814) q[0];
rz(3.0389298) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(1.7662802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1840399) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(-2.1368105) q[1];
rz(-pi) q[2];
rz(0.41820742) q[3];
sx q[3];
rz(-1.4718664) q[3];
sx q[3];
rz(-1.4560771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2404279) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(3.124776) q[0];
rz(3.1230714) q[1];
sx q[1];
rz(-1.3341981) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803857) q[0];
sx q[0];
rz(-1.2347504) q[0];
sx q[0];
rz(1.8340322) q[0];
rz(-pi) q[1];
rz(0.2741371) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(2.0704839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2494431) q[1];
sx q[1];
rz(-1.8037233) q[1];
sx q[1];
rz(-1.1942785) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7917463) q[3];
sx q[3];
rz(-1.5077935) q[3];
sx q[3];
rz(2.8642879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(2.418628) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13847362) q[0];
sx q[0];
rz(-1.4644633) q[0];
sx q[0];
rz(-1.7220201) q[0];
rz(-pi) q[1];
rz(0.16913551) q[2];
sx q[2];
rz(-1.9431056) q[2];
sx q[2];
rz(1.7793836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1180229) q[1];
sx q[1];
rz(-1.6260864) q[1];
sx q[1];
rz(2.0748595) q[1];
rz(-pi) q[2];
x q[2];
rz(2.74182) q[3];
sx q[3];
rz(-2.3120566) q[3];
sx q[3];
rz(1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4621949) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(0.47808403) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(0.59393926) q[3];
sx q[3];
rz(-2.754302) q[3];
sx q[3];
rz(2.1828628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
