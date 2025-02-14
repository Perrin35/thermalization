OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0107083) q[0];
sx q[0];
rz(-2.4680128) q[0];
sx q[0];
rz(-2.3186865) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17774489) q[0];
sx q[0];
rz(-1.2262934) q[0];
sx q[0];
rz(2.3314407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0124341) q[2];
sx q[2];
rz(-1.2744546) q[2];
sx q[2];
rz(1.3989965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.365578) q[1];
sx q[1];
rz(-1.6322337) q[1];
sx q[1];
rz(-0.12428026) q[1];
rz(-3.0700686) q[3];
sx q[3];
rz(-1.2010964) q[3];
sx q[3];
rz(2.2184851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9660008) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(-1.9123745) q[2];
rz(1.2980596) q[3];
sx q[3];
rz(-1.3763873) q[3];
sx q[3];
rz(-1.957533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-0.17756473) q[0];
sx q[0];
rz(-2.8819045) q[0];
sx q[0];
rz(0.92726707) q[0];
rz(-0.19993965) q[1];
sx q[1];
rz(-1.6688469) q[1];
sx q[1];
rz(-2.1520069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703116) q[0];
sx q[0];
rz(-1.6785445) q[0];
sx q[0];
rz(-1.2077483) q[0];
rz(-pi) q[1];
rz(-1.6140429) q[2];
sx q[2];
rz(-2.0048365) q[2];
sx q[2];
rz(-2.1021646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4633219) q[1];
sx q[1];
rz(-0.090677977) q[1];
sx q[1];
rz(1.00738) q[1];
x q[2];
rz(2.5603676) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(0.083409781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9280615) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(-2.7878063) q[2];
rz(-0.95156041) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(0.32683867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3554409) q[0];
sx q[0];
rz(-2.1936301) q[0];
sx q[0];
rz(1.0107262) q[0];
rz(3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-0.40930632) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362013) q[0];
sx q[0];
rz(-2.276148) q[0];
sx q[0];
rz(-0.58569293) q[0];
x q[1];
rz(-1.1597761) q[2];
sx q[2];
rz(-1.4189229) q[2];
sx q[2];
rz(1.077026) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61534897) q[1];
sx q[1];
rz(-2.1110635) q[1];
sx q[1];
rz(-1.8002285) q[1];
rz(-1.5915967) q[3];
sx q[3];
rz(-2.0796418) q[3];
sx q[3];
rz(1.5326981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.056444082) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(0.15963456) q[2];
rz(1.2498445) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(-1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98642629) q[0];
sx q[0];
rz(-0.43743375) q[0];
sx q[0];
rz(-1.0945818) q[0];
rz(1.1711586) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(-1.0135244) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60166925) q[0];
sx q[0];
rz(-3.0832096) q[0];
sx q[0];
rz(-2.4235382) q[0];
rz(-1.0763747) q[2];
sx q[2];
rz(-1.5740047) q[2];
sx q[2];
rz(2.6319844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2029548) q[1];
sx q[1];
rz(-1.5302916) q[1];
sx q[1];
rz(0.26085965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7988241) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(2.5955615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0052428) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(1.8278149) q[2];
rz(-0.93787307) q[3];
sx q[3];
rz(-1.2910941) q[3];
sx q[3];
rz(-3.042799) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(2.0966356) q[0];
rz(0.16432556) q[1];
sx q[1];
rz(-1.798809) q[1];
sx q[1];
rz(-2.9716861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2192046) q[0];
sx q[0];
rz(-1.4996075) q[0];
sx q[0];
rz(1.7573331) q[0];
rz(-pi) q[1];
rz(-0.33268945) q[2];
sx q[2];
rz(-1.2342808) q[2];
sx q[2];
rz(2.3115186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8039717) q[1];
sx q[1];
rz(-0.70120431) q[1];
sx q[1];
rz(-1.1887202) q[1];
rz(-pi) q[2];
rz(-1.4395797) q[3];
sx q[3];
rz(-2.2926712) q[3];
sx q[3];
rz(-0.038824507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56875151) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(-1.0315726) q[2];
rz(-2.0555563) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.731855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373435) q[0];
sx q[0];
rz(-1.046109) q[0];
sx q[0];
rz(1.401249) q[0];
rz(-0.42964545) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(-1.8362129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25241643) q[0];
sx q[0];
rz(-0.85491291) q[0];
sx q[0];
rz(1.9837911) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1370509) q[2];
sx q[2];
rz(-1.1318739) q[2];
sx q[2];
rz(0.86290854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67638799) q[1];
sx q[1];
rz(-2.891245) q[1];
sx q[1];
rz(-3.0167104) q[1];
rz(-pi) q[2];
rz(-1.3280198) q[3];
sx q[3];
rz(-1.4221898) q[3];
sx q[3];
rz(-2.7000303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0171011) q[2];
sx q[2];
rz(-1.93511) q[2];
sx q[2];
rz(-0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5519578) q[3];
sx q[3];
rz(2.8563833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8816836) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(0.024070865) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(-0.37758652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51272362) q[0];
sx q[0];
rz(-2.2170904) q[0];
sx q[0];
rz(-1.76684) q[0];
x q[1];
rz(1.3578916) q[2];
sx q[2];
rz(-1.7434374) q[2];
sx q[2];
rz(1.5162692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3096034) q[1];
sx q[1];
rz(-1.9506694) q[1];
sx q[1];
rz(-1.1516476) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33414109) q[3];
sx q[3];
rz(-1.982882) q[3];
sx q[3];
rz(-2.0711632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5923803) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(0.40360061) q[2];
rz(2.9523383) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(0.45690593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6117578) q[0];
sx q[0];
rz(-0.54946041) q[0];
sx q[0];
rz(2.8305565) q[0];
rz(0.95651904) q[1];
sx q[1];
rz(-1.6638959) q[1];
sx q[1];
rz(0.32435736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7081147) q[0];
sx q[0];
rz(-1.38033) q[0];
sx q[0];
rz(-2.6300927) q[0];
x q[1];
rz(1.5115405) q[2];
sx q[2];
rz(-1.0612592) q[2];
sx q[2];
rz(-1.6409724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2203103) q[1];
sx q[1];
rz(-1.1571572) q[1];
sx q[1];
rz(2.2181554) q[1];
rz(-pi) q[2];
rz(-2.0288386) q[3];
sx q[3];
rz(-2.4107217) q[3];
sx q[3];
rz(-2.3422086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73584622) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(-0.22845593) q[2];
rz(0.53260803) q[3];
sx q[3];
rz(-0.76627982) q[3];
sx q[3];
rz(-2.2920091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547895) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(1.2702031) q[0];
rz(-0.58865976) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(0.38280815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7080806) q[0];
sx q[0];
rz(-1.5956586) q[0];
sx q[0];
rz(-1.8000425) q[0];
rz(-pi) q[1];
rz(-0.43395502) q[2];
sx q[2];
rz(-1.343691) q[2];
sx q[2];
rz(-1.3004829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36675699) q[1];
sx q[1];
rz(-2.2765571) q[1];
sx q[1];
rz(-1.4053759) q[1];
rz(-pi) q[2];
rz(2.8715897) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(-1.6635513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(-2.7413979) q[2];
rz(-2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(-2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8409214) q[0];
sx q[0];
rz(-1.3341757) q[0];
sx q[0];
rz(1.0155431) q[0];
rz(-0.66889846) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(1.5060172) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5251314) q[0];
sx q[0];
rz(-0.86989738) q[0];
sx q[0];
rz(2.4674795) q[0];
rz(-pi) q[1];
rz(2.593562) q[2];
sx q[2];
rz(-1.9491674) q[2];
sx q[2];
rz(0.038383287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4989711) q[1];
sx q[1];
rz(-1.6097415) q[1];
sx q[1];
rz(-1.0907111) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90088441) q[3];
sx q[3];
rz(-2.3436119) q[3];
sx q[3];
rz(-0.98679286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(1.6157185) q[2];
rz(-2.8499991) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780554) q[0];
sx q[0];
rz(-0.96374496) q[0];
sx q[0];
rz(-2.5441334) q[0];
rz(-2.8529104) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(2.504442) q[2];
sx q[2];
rz(-0.49497866) q[2];
sx q[2];
rz(-0.53866932) q[2];
rz(-0.47824974) q[3];
sx q[3];
rz(-0.50898715) q[3];
sx q[3];
rz(-0.57686808) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
