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
rz(-1.7614814) q[0];
sx q[0];
rz(-1.5505646) q[0];
sx q[0];
rz(0.38183364) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782717) q[0];
sx q[0];
rz(-2.1594783) q[0];
sx q[0];
rz(0.6529385) q[0];
rz(-pi) q[1];
rz(-1.5906232) q[2];
sx q[2];
rz(-1.9409436) q[2];
sx q[2];
rz(-1.4408979) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4630064) q[1];
sx q[1];
rz(-1.9941115) q[1];
sx q[1];
rz(1.4284575) q[1];
rz(1.0409058) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(-1.4640946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(-1.5748242) q[2];
rz(-1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(-0.69275698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67398706) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(3.0740956) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(0.31879058) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2683811) q[0];
sx q[0];
rz(-2.2751774) q[0];
sx q[0];
rz(-1.1581139) q[0];
rz(-pi) q[1];
rz(2.6574357) q[2];
sx q[2];
rz(-0.47704298) q[2];
sx q[2];
rz(2.8271528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31736703) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(2.1721858) q[1];
rz(-pi) q[2];
rz(3.0158668) q[3];
sx q[3];
rz(-1.8549524) q[3];
sx q[3];
rz(-1.3181669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(-0.47743615) q[2];
rz(-1.3580953) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(0.40621743) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3258404) q[0];
sx q[0];
rz(-1.6597972) q[0];
sx q[0];
rz(-2.1584216) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94495541) q[2];
sx q[2];
rz(-1.8650818) q[2];
sx q[2];
rz(-1.6950032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.660608) q[1];
sx q[1];
rz(-2.0913731) q[1];
sx q[1];
rz(2.1064227) q[1];
rz(-pi) q[2];
rz(2.0800679) q[3];
sx q[3];
rz(-1.2767999) q[3];
sx q[3];
rz(-1.7727838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4262126) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(1.0769843) q[2];
rz(1.2601323) q[3];
sx q[3];
rz(-2.1144805) q[3];
sx q[3];
rz(0.25585678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28427163) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(-2.381109) q[0];
rz(-2.7898232) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(-0.15028353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6929157) q[0];
sx q[0];
rz(-1.5741072) q[0];
sx q[0];
rz(3.1396554) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1614805) q[2];
sx q[2];
rz(-2.7285353) q[2];
sx q[2];
rz(1.6020136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3327738) q[1];
sx q[1];
rz(-1.2398232) q[1];
sx q[1];
rz(0.95167758) q[1];
rz(-pi) q[2];
rz(1.3129381) q[3];
sx q[3];
rz(-1.1262731) q[3];
sx q[3];
rz(0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.030423) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(2.5566768) q[2];
rz(3.0758744) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(1.0544624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86972648) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(-0.72750339) q[0];
rz(-1.2959405) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(-2.4268699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.489577) q[0];
sx q[0];
rz(-1.4997109) q[0];
sx q[0];
rz(-2.9201304) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7834846) q[2];
sx q[2];
rz(-0.3198959) q[2];
sx q[2];
rz(1.7826155) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2421869) q[1];
sx q[1];
rz(-1.9493628) q[1];
sx q[1];
rz(-0.19902163) q[1];
x q[2];
rz(-0.994579) q[3];
sx q[3];
rz(-1.2354697) q[3];
sx q[3];
rz(-1.9063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99594816) q[2];
sx q[2];
rz(-1.5278634) q[2];
sx q[2];
rz(-1.9963416) q[2];
rz(2.9389985) q[3];
sx q[3];
rz(-0.069772094) q[3];
sx q[3];
rz(0.32874671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458493) q[0];
sx q[0];
rz(-2.8491617) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(2.5766418) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(-1.8126743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5707626) q[0];
sx q[0];
rz(-2.0432297) q[0];
sx q[0];
rz(-2.5357312) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6199048) q[2];
sx q[2];
rz(-1.7108166) q[2];
sx q[2];
rz(-2.5852709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5766474) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(0.74585657) q[1];
rz(-2.9181913) q[3];
sx q[3];
rz(-1.3376682) q[3];
sx q[3];
rz(0.65321556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38587511) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(2.7937549) q[2];
rz(-2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(-0.95112479) q[0];
rz(-2.8490207) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(0.94400418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2068588) q[0];
sx q[0];
rz(-0.30843302) q[0];
sx q[0];
rz(-1.5200204) q[0];
x q[1];
rz(-2.6636276) q[2];
sx q[2];
rz(-1.1349808) q[2];
sx q[2];
rz(-2.8675241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4237562) q[1];
sx q[1];
rz(-1.5482117) q[1];
sx q[1];
rz(-1.3976239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1066426) q[3];
sx q[3];
rz(-0.45848819) q[3];
sx q[3];
rz(1.637984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49147478) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(0.024638351) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(-2.1520481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330419) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(-3.0035875) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-2.2053351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662921) q[0];
sx q[0];
rz(-1.202824) q[0];
sx q[0];
rz(1.8938246) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8291275) q[2];
sx q[2];
rz(-1.1489604) q[2];
sx q[2];
rz(1.2527996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1106134) q[1];
sx q[1];
rz(-0.18337164) q[1];
sx q[1];
rz(-1.490209) q[1];
rz(-pi) q[2];
rz(0.68558399) q[3];
sx q[3];
rz(-1.778025) q[3];
sx q[3];
rz(-0.7570467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.998698) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(-0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1634624) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(1.7568024) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(-0.91420954) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5152138) q[0];
sx q[0];
rz(-2.6914586) q[0];
sx q[0];
rz(-1.2594957) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73458521) q[2];
sx q[2];
rz(-1.6050287) q[2];
sx q[2];
rz(1.5058668) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2466941) q[1];
sx q[1];
rz(-1.7082038) q[1];
sx q[1];
rz(-0.64260428) q[1];
rz(-pi) q[2];
rz(-2.8960885) q[3];
sx q[3];
rz(-2.5872018) q[3];
sx q[3];
rz(-1.7911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(1.6287104) q[2];
rz(-0.58879876) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-2.5741757) q[0];
sx q[0];
rz(-0.2844511) q[0];
rz(-0.65471634) q[1];
sx q[1];
rz(-1.3755211) q[1];
sx q[1];
rz(0.42025748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81440825) q[0];
sx q[0];
rz(-1.1564009) q[0];
sx q[0];
rz(-0.23420686) q[0];
x q[1];
rz(-2.0345694) q[2];
sx q[2];
rz(-2.3056185) q[2];
sx q[2];
rz(-2.7157264) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79607591) q[1];
sx q[1];
rz(-2.9931167) q[1];
sx q[1];
rz(-2.2125437) q[1];
rz(-pi) q[2];
rz(-2.5758366) q[3];
sx q[3];
rz(-1.2671736) q[3];
sx q[3];
rz(-2.2591801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(-2.1238756) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(-0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9781072) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(0.12374395) q[1];
sx q[1];
rz(-2.1950304) q[1];
sx q[1];
rz(-1.310941) q[1];
rz(1.704557) q[2];
sx q[2];
rz(-1.8885713) q[2];
sx q[2];
rz(1.0071913) q[2];
rz(-0.10931482) q[3];
sx q[3];
rz(-0.31382618) q[3];
sx q[3];
rz(-1.8574497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
