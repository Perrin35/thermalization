OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(8.6241047) q[0];
sx q[0];
rz(9.5700349) q[0];
rz(-0.89291209) q[1];
sx q[1];
rz(3.555759) q[1];
sx q[1];
rz(10.531737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205638) q[0];
sx q[0];
rz(-2.1269089) q[0];
sx q[0];
rz(-0.70379852) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28319226) q[2];
sx q[2];
rz(-1.257466) q[2];
sx q[2];
rz(1.6746132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2295215) q[1];
sx q[1];
rz(-2.1451277) q[1];
sx q[1];
rz(-0.34194754) q[1];
rz(-2.0756196) q[3];
sx q[3];
rz(-2.7144066) q[3];
sx q[3];
rz(1.2963201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5408111) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-3.0493128) q[2];
rz(-1.4394834) q[3];
sx q[3];
rz(-1.9974134) q[3];
sx q[3];
rz(0.93851844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.1630031) q[0];
sx q[0];
rz(2.5376885) q[0];
rz(1.6101135) q[1];
sx q[1];
rz(-1.3319301) q[1];
sx q[1];
rz(-1.5276705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85717843) q[0];
sx q[0];
rz(-2.1499794) q[0];
sx q[0];
rz(0.37838899) q[0];
x q[1];
rz(0.98601922) q[2];
sx q[2];
rz(-0.9169609) q[2];
sx q[2];
rz(-0.80160917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70992596) q[1];
sx q[1];
rz(-1.559169) q[1];
sx q[1];
rz(-0.27375582) q[1];
x q[2];
rz(-0.098193125) q[3];
sx q[3];
rz(-1.5714688) q[3];
sx q[3];
rz(2.0405586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4332726) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(1.0478919) q[2];
rz(-0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3480551) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(1.0070356) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(-0.67063355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8438727) q[0];
sx q[0];
rz(-0.62652912) q[0];
sx q[0];
rz(-1.7202176) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3449981) q[2];
sx q[2];
rz(-2.0907986) q[2];
sx q[2];
rz(0.3546302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.933316) q[1];
sx q[1];
rz(-1.2348935) q[1];
sx q[1];
rz(-0.64888727) q[1];
rz(1.8995883) q[3];
sx q[3];
rz(-0.34600779) q[3];
sx q[3];
rz(0.11869988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.85155073) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(-2.2500989) q[2];
rz(-1.1024891) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(1.4055143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3201228) q[0];
sx q[0];
rz(-1.3865043) q[0];
sx q[0];
rz(-2.084305) q[0];
rz(3.140246) q[1];
sx q[1];
rz(-2.3206382) q[1];
sx q[1];
rz(-2.3473306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43633902) q[0];
sx q[0];
rz(-2.2326075) q[0];
sx q[0];
rz(-1.5474942) q[0];
x q[1];
rz(-2.4636872) q[2];
sx q[2];
rz(-2.6982582) q[2];
sx q[2];
rz(1.5405129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92137488) q[1];
sx q[1];
rz(-1.1168606) q[1];
sx q[1];
rz(-2.3364319) q[1];
x q[2];
rz(-2.7806588) q[3];
sx q[3];
rz(-2.8246701) q[3];
sx q[3];
rz(1.288687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0346251) q[2];
sx q[2];
rz(-2.5204973) q[2];
sx q[2];
rz(1.0221488) q[2];
rz(-1.243783) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(-1.3589842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-2.6720402) q[0];
sx q[0];
rz(-1.38009) q[0];
sx q[0];
rz(0.45355466) q[0];
rz(-2.1039311) q[1];
sx q[1];
rz(-1.1353759) q[1];
sx q[1];
rz(-0.79197788) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29402367) q[0];
sx q[0];
rz(-1.5269465) q[0];
sx q[0];
rz(1.8490318) q[0];
rz(-pi) q[1];
rz(1.0019962) q[2];
sx q[2];
rz(-2.482907) q[2];
sx q[2];
rz(1.9994761) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5752069) q[1];
sx q[1];
rz(-0.35950152) q[1];
sx q[1];
rz(2.8043037) q[1];
rz(-pi) q[2];
rz(1.7030992) q[3];
sx q[3];
rz(-3.0022394) q[3];
sx q[3];
rz(2.9840368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8301293) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(0.30544454) q[2];
rz(-2.9597802) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.55564725) q[0];
sx q[0];
rz(-2.3190627) q[0];
sx q[0];
rz(-0.12538759) q[0];
rz(-1.5665945) q[1];
sx q[1];
rz(-1.6927203) q[1];
sx q[1];
rz(-3.1094508) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956995) q[0];
sx q[0];
rz(-2.3902378) q[0];
sx q[0];
rz(0.32352792) q[0];
rz(-pi) q[1];
rz(1.6400385) q[2];
sx q[2];
rz(-2.4089795) q[2];
sx q[2];
rz(-2.46012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5405724) q[1];
sx q[1];
rz(-2.2443612) q[1];
sx q[1];
rz(2.5190398) q[1];
rz(-pi) q[2];
rz(-0.7456268) q[3];
sx q[3];
rz(-1.485482) q[3];
sx q[3];
rz(-1.2554393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0142168) q[2];
sx q[2];
rz(-1.9275503) q[2];
sx q[2];
rz(-2.0188913) q[2];
rz(-3.0681916) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(0.61298031) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70925322) q[0];
sx q[0];
rz(-1.657635) q[0];
sx q[0];
rz(-0.51026979) q[0];
rz(0.36422745) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(1.4998923) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071935805) q[0];
sx q[0];
rz(-2.6652415) q[0];
sx q[0];
rz(1.8029965) q[0];
rz(-pi) q[1];
rz(-1.4128311) q[2];
sx q[2];
rz(-1.4539495) q[2];
sx q[2];
rz(-0.37010461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81699521) q[1];
sx q[1];
rz(-2.5909068) q[1];
sx q[1];
rz(2.703642) q[1];
rz(-pi) q[2];
rz(-0.73697258) q[3];
sx q[3];
rz(-1.342257) q[3];
sx q[3];
rz(1.6933954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69786543) q[2];
sx q[2];
rz(-1.5550193) q[2];
sx q[2];
rz(-0.86722428) q[2];
rz(-1.5813658) q[3];
sx q[3];
rz(-0.19872228) q[3];
sx q[3];
rz(1.3377415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7241868) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(0.41626406) q[0];
rz(-1.9789713) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(2.3847041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1838186) q[0];
sx q[0];
rz(-1.8584434) q[0];
sx q[0];
rz(-1.1270866) q[0];
rz(-pi) q[1];
rz(-0.9932809) q[2];
sx q[2];
rz(-1.9406291) q[2];
sx q[2];
rz(0.12060697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97786602) q[1];
sx q[1];
rz(-1.1440579) q[1];
sx q[1];
rz(0.20789288) q[1];
rz(-pi) q[2];
rz(-1.105433) q[3];
sx q[3];
rz(-0.98794395) q[3];
sx q[3];
rz(0.40377709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4550712) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(-1.7835468) q[2];
rz(-0.81651917) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(2.9787298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002976) q[0];
sx q[0];
rz(-1.4485899) q[0];
sx q[0];
rz(-1.7342389) q[0];
rz(0.74527144) q[1];
sx q[1];
rz(-0.73423568) q[1];
sx q[1];
rz(1.5596681) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1187706) q[0];
sx q[0];
rz(-1.2948991) q[0];
sx q[0];
rz(-2.2488382) q[0];
x q[1];
rz(-3.09868) q[2];
sx q[2];
rz(-2.4933698) q[2];
sx q[2];
rz(-0.15951482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.078439039) q[1];
sx q[1];
rz(-2.6780824) q[1];
sx q[1];
rz(1.2786675) q[1];
rz(-pi) q[2];
rz(-1.6635062) q[3];
sx q[3];
rz(-2.5639236) q[3];
sx q[3];
rz(0.81837624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72426307) q[2];
sx q[2];
rz(-1.4549007) q[2];
sx q[2];
rz(1.1154741) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-0.0089664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(2.3790835) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-1.0647048) q[1];
sx q[1];
rz(-1.9047838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002703) q[0];
sx q[0];
rz(-1.5093439) q[0];
sx q[0];
rz(0.49646722) q[0];
rz(-pi) q[1];
rz(-2.2549596) q[2];
sx q[2];
rz(-1.1012474) q[2];
sx q[2];
rz(-2.4327715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0227388) q[1];
sx q[1];
rz(-0.035422649) q[1];
sx q[1];
rz(-0.78411786) q[1];
rz(-pi) q[2];
rz(-1.4051626) q[3];
sx q[3];
rz(-2.7319623) q[3];
sx q[3];
rz(-1.7543242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2463871) q[2];
sx q[2];
rz(-3.0463986) q[2];
sx q[2];
rz(3.134356) q[2];
rz(0.17035189) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(-0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728747) q[0];
sx q[0];
rz(-1.3997411) q[0];
sx q[0];
rz(1.1135143) q[0];
rz(3.0524104) q[1];
sx q[1];
rz(-1.6249648) q[1];
sx q[1];
rz(-1.9393495) q[1];
rz(-0.15684814) q[2];
sx q[2];
rz(-1.7050171) q[2];
sx q[2];
rz(-0.27324054) q[2];
rz(-0.54696541) q[3];
sx q[3];
rz(-0.34815683) q[3];
sx q[3];
rz(1.9937594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
