OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.469874) q[0];
sx q[0];
rz(-1.9865541) q[0];
sx q[0];
rz(1.7116829) q[0];
rz(0.45541304) q[1];
sx q[1];
rz(-0.59133426) q[1];
sx q[1];
rz(2.0145745) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6756962) q[0];
sx q[0];
rz(-1.5592335) q[0];
sx q[0];
rz(1.3402864) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65140407) q[2];
sx q[2];
rz(-0.42170516) q[2];
sx q[2];
rz(2.5159474) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.110636) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(0.33302506) q[1];
rz(-pi) q[2];
rz(-1.5886878) q[3];
sx q[3];
rz(-1.3800637) q[3];
sx q[3];
rz(2.408556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0226125) q[2];
sx q[2];
rz(-2.4369414) q[2];
sx q[2];
rz(-1.653778) q[2];
rz(0.37114272) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356075) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(2.4775179) q[0];
rz(-2.5750419) q[1];
sx q[1];
rz(-1.605875) q[1];
sx q[1];
rz(0.18633349) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5390426) q[0];
sx q[0];
rz(-1.3485753) q[0];
sx q[0];
rz(1.0265539) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9695942) q[2];
sx q[2];
rz(-2.6896713) q[2];
sx q[2];
rz(-2.4659803) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1809606) q[1];
sx q[1];
rz(-0.69968984) q[1];
sx q[1];
rz(-1.3153428) q[1];
x q[2];
rz(1.6870895) q[3];
sx q[3];
rz(-0.52403677) q[3];
sx q[3];
rz(-0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(1.5796278) q[2];
rz(-1.1473848) q[3];
sx q[3];
rz(-0.29427823) q[3];
sx q[3];
rz(1.7779721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550302) q[0];
sx q[0];
rz(-1.6494305) q[0];
sx q[0];
rz(2.8360039) q[0];
rz(-1.8606404) q[1];
sx q[1];
rz(-0.31550229) q[1];
sx q[1];
rz(-1.5234647) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9010503) q[0];
sx q[0];
rz(-0.43473703) q[0];
sx q[0];
rz(2.5655599) q[0];
rz(-pi) q[1];
rz(0.32714897) q[2];
sx q[2];
rz(-0.047657813) q[2];
sx q[2];
rz(0.41027712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2669526) q[1];
sx q[1];
rz(-0.8921418) q[1];
sx q[1];
rz(2.0372422) q[1];
x q[2];
rz(1.9033857) q[3];
sx q[3];
rz(-2.5032515) q[3];
sx q[3];
rz(2.7901543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65172255) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(-2.9215802) q[2];
rz(3.0618727) q[3];
sx q[3];
rz(-2.1646175) q[3];
sx q[3];
rz(2.2102833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109969) q[0];
sx q[0];
rz(-1.7544704) q[0];
sx q[0];
rz(3.1331449) q[0];
rz(1.1620109) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(2.0268424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.312437) q[0];
sx q[0];
rz(-0.36530802) q[0];
sx q[0];
rz(2.1213953) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1276541) q[2];
sx q[2];
rz(-1.077855) q[2];
sx q[2];
rz(-1.8007474) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1871918) q[1];
sx q[1];
rz(-1.2174509) q[1];
sx q[1];
rz(0.29275972) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1590121) q[3];
sx q[3];
rz(-1.1819541) q[3];
sx q[3];
rz(1.6568615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70725012) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(0.84558359) q[2];
rz(-0.25104684) q[3];
sx q[3];
rz(-1.2716525) q[3];
sx q[3];
rz(-2.3959851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71395981) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(1.0300256) q[0];
rz(-2.5469942) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(1.7880218) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1270238) q[0];
sx q[0];
rz(-1.5935803) q[0];
sx q[0];
rz(1.5915074) q[0];
rz(0.51026235) q[2];
sx q[2];
rz(-2.7760421) q[2];
sx q[2];
rz(0.22399677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7673257) q[1];
sx q[1];
rz(-2.0637636) q[1];
sx q[1];
rz(-0.79549353) q[1];
x q[2];
rz(1.0566125) q[3];
sx q[3];
rz(-2.100889) q[3];
sx q[3];
rz(-0.3584273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2305962) q[2];
sx q[2];
rz(-1.6350919) q[2];
sx q[2];
rz(-0.050749151) q[2];
rz(1.1132318) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(0.39065233) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0942866) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(-0.37288368) q[0];
rz(-1.3039543) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(2.1162927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4639125) q[0];
sx q[0];
rz(-0.45651528) q[0];
sx q[0];
rz(2.9677261) q[0];
rz(-0.16496678) q[2];
sx q[2];
rz(-1.2447164) q[2];
sx q[2];
rz(0.13298377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2220425) q[1];
sx q[1];
rz(-1.371438) q[1];
sx q[1];
rz(0.87032236) q[1];
rz(-1.9944836) q[3];
sx q[3];
rz(-1.7762587) q[3];
sx q[3];
rz(0.68483554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9103526) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.4300038) q[2];
rz(-1.5405687) q[3];
sx q[3];
rz(-2.4547596) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533399) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(0.43149313) q[0];
rz(1.2639812) q[1];
sx q[1];
rz(-1.6004205) q[1];
sx q[1];
rz(-0.65013179) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6912731) q[0];
sx q[0];
rz(-0.29940816) q[0];
sx q[0];
rz(2.8595238) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88754295) q[2];
sx q[2];
rz(-1.4897963) q[2];
sx q[2];
rz(1.1095059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7909779) q[1];
sx q[1];
rz(-0.87360604) q[1];
sx q[1];
rz(-0.15048262) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1339478) q[3];
sx q[3];
rz(-1.1129324) q[3];
sx q[3];
rz(-1.0308105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55107895) q[2];
sx q[2];
rz(-1.3912018) q[2];
sx q[2];
rz(3.1371269) q[2];
rz(-3.1320069) q[3];
sx q[3];
rz(-1.8066581) q[3];
sx q[3];
rz(2.7981304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.9850605) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(0.16047934) q[0];
rz(-1.7566682) q[1];
sx q[1];
rz(-0.79912186) q[1];
sx q[1];
rz(2.8327732) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35058103) q[0];
sx q[0];
rz(-1.6948943) q[0];
sx q[0];
rz(-1.3826293) q[0];
x q[1];
rz(0.71134348) q[2];
sx q[2];
rz(-1.8677731) q[2];
sx q[2];
rz(-2.1893196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2330197) q[1];
sx q[1];
rz(-1.1814139) q[1];
sx q[1];
rz(-1.1638457) q[1];
x q[2];
rz(-0.35328226) q[3];
sx q[3];
rz(-1.8646514) q[3];
sx q[3];
rz(-2.053956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7979692) q[2];
sx q[2];
rz(-2.5324731) q[2];
sx q[2];
rz(-1.5416175) q[2];
rz(-2.4566417) q[3];
sx q[3];
rz(-1.5870321) q[3];
sx q[3];
rz(-1.6182914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63055444) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(0.55139971) q[0];
rz(-0.4772056) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(-2.3068857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4906368) q[0];
sx q[0];
rz(-1.0792562) q[0];
sx q[0];
rz(-2.0975344) q[0];
rz(-0.70087011) q[2];
sx q[2];
rz(-1.5972023) q[2];
sx q[2];
rz(0.7757265) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0352201) q[1];
sx q[1];
rz(-2.1823931) q[1];
sx q[1];
rz(-0.04208438) q[1];
rz(-pi) q[2];
rz(2.5284994) q[3];
sx q[3];
rz(-1.0524155) q[3];
sx q[3];
rz(-1.6676774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4386374) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(2.9023602) q[2];
rz(-0.3244102) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(1.602406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5156373) q[0];
sx q[0];
rz(-3.0037168) q[0];
sx q[0];
rz(2.4249518) q[0];
rz(0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(0.31002054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6577756) q[0];
sx q[0];
rz(-0.80379009) q[0];
sx q[0];
rz(2.0652384) q[0];
rz(-2.7298195) q[2];
sx q[2];
rz(-1.207282) q[2];
sx q[2];
rz(-2.1902254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2523732) q[1];
sx q[1];
rz(-2.3205726) q[1];
sx q[1];
rz(2.0682425) q[1];
rz(-pi) q[2];
rz(2.4118092) q[3];
sx q[3];
rz(-0.68396362) q[3];
sx q[3];
rz(1.3866977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50905716) q[2];
sx q[2];
rz(-0.57764235) q[2];
sx q[2];
rz(1.4198111) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-0.71075478) q[3];
sx q[3];
rz(2.3783309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2390908) q[0];
sx q[0];
rz(-2.1949174) q[0];
sx q[0];
rz(1.3876023) q[0];
rz(-1.6613962) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(-2.8724332) q[2];
sx q[2];
rz(-2.8648389) q[2];
sx q[2];
rz(1.7493389) q[2];
rz(2.1275413) q[3];
sx q[3];
rz(-0.76273668) q[3];
sx q[3];
rz(-0.098165011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
