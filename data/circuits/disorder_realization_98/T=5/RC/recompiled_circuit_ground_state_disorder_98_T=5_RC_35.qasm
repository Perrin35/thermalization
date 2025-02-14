OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7336361) q[0];
sx q[0];
rz(-0.1853369) q[0];
sx q[0];
rz(1.4987401) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(2.763881) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2687338) q[0];
sx q[0];
rz(-1.3493363) q[0];
sx q[0];
rz(-2.8941157) q[0];
x q[1];
rz(2.0091363) q[2];
sx q[2];
rz(-1.1852263) q[2];
sx q[2];
rz(1.350268) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7491748) q[1];
sx q[1];
rz(-0.91411829) q[1];
sx q[1];
rz(1.6139469) q[1];
rz(-pi) q[2];
rz(1.6781647) q[3];
sx q[3];
rz(-0.26623785) q[3];
sx q[3];
rz(-0.050198089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3898042) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(2.8765615) q[2];
rz(0.42936471) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.4877743) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(-1.8408467) q[0];
rz(3.0744413) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(-0.86004177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9275151) q[0];
sx q[0];
rz(-0.96827114) q[0];
sx q[0];
rz(-2.430315) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2566893) q[2];
sx q[2];
rz(-0.98457789) q[2];
sx q[2];
rz(-0.84916678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1971284) q[1];
sx q[1];
rz(-1.010506) q[1];
sx q[1];
rz(-0.48546882) q[1];
rz(1.0299171) q[3];
sx q[3];
rz(-2.7949998) q[3];
sx q[3];
rz(2.5427713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63090515) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(0.56646937) q[2];
rz(2.412292) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(-0.20310371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.94839621) q[0];
sx q[0];
rz(-2.0553135) q[0];
sx q[0];
rz(-2.7753944) q[0];
rz(1.5557479) q[1];
sx q[1];
rz(-2.0987174) q[1];
sx q[1];
rz(0.050447024) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4719066) q[0];
sx q[0];
rz(-1.4916873) q[0];
sx q[0];
rz(-3.0997653) q[0];
rz(1.0773727) q[2];
sx q[2];
rz(-0.9443379) q[2];
sx q[2];
rz(2.3979417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55712426) q[1];
sx q[1];
rz(-1.4038175) q[1];
sx q[1];
rz(-2.9191769) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8837711) q[3];
sx q[3];
rz(-0.44383263) q[3];
sx q[3];
rz(2.7601506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(-1.3107497) q[2];
rz(-1.358076) q[3];
sx q[3];
rz(-2.570593) q[3];
sx q[3];
rz(-1.8324435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27291372) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(1.7568461) q[0];
rz(-1.9028496) q[1];
sx q[1];
rz(-1.2307931) q[1];
sx q[1];
rz(1.1741656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096701972) q[0];
sx q[0];
rz(-2.6706123) q[0];
sx q[0];
rz(1.1440008) q[0];
rz(0.73100369) q[2];
sx q[2];
rz(-0.36492294) q[2];
sx q[2];
rz(0.66563767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3191022) q[1];
sx q[1];
rz(-0.61885364) q[1];
sx q[1];
rz(-0.98554218) q[1];
x q[2];
rz(-2.2635095) q[3];
sx q[3];
rz(-1.2125848) q[3];
sx q[3];
rz(-2.3830551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43526402) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(3.1026133) q[2];
rz(0.6428166) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(2.5207998) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3103127) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(0.020462791) q[0];
rz(2.9488355) q[1];
sx q[1];
rz(-2.5242476) q[1];
sx q[1];
rz(-1.1654759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63686164) q[0];
sx q[0];
rz(-0.76380542) q[0];
sx q[0];
rz(-2.5434407) q[0];
rz(-pi) q[1];
rz(-2.5978659) q[2];
sx q[2];
rz(-2.3909878) q[2];
sx q[2];
rz(0.0074530938) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36042696) q[1];
sx q[1];
rz(-1.5910223) q[1];
sx q[1];
rz(1.0670877) q[1];
x q[2];
rz(-1.9475157) q[3];
sx q[3];
rz(-2.93749) q[3];
sx q[3];
rz(-0.78777504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8274902) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(2.6143383) q[2];
rz(-0.54316795) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813893) q[0];
sx q[0];
rz(-0.15601604) q[0];
sx q[0];
rz(-2.7348837) q[0];
rz(-1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(-1.0103753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8445379) q[0];
sx q[0];
rz(-2.5803891) q[0];
sx q[0];
rz(-1.7991174) q[0];
x q[1];
rz(0.40596227) q[2];
sx q[2];
rz(-0.34892198) q[2];
sx q[2];
rz(2.5418848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53394371) q[1];
sx q[1];
rz(-1.2022809) q[1];
sx q[1];
rz(1.7979969) q[1];
x q[2];
rz(-2.3361037) q[3];
sx q[3];
rz(-1.3873161) q[3];
sx q[3];
rz(0.047540548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0747718) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(2.1066966) q[2];
rz(-1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9676301) q[0];
sx q[0];
rz(-1.5246464) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(-0.11218849) q[1];
sx q[1];
rz(-1.2022377) q[1];
sx q[1];
rz(-0.97253886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93720651) q[0];
sx q[0];
rz(-2.4582259) q[0];
sx q[0];
rz(-2.1478189) q[0];
rz(0.99268408) q[2];
sx q[2];
rz(-0.51380605) q[2];
sx q[2];
rz(-2.2058599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9303846) q[1];
sx q[1];
rz(-0.58961419) q[1];
sx q[1];
rz(0.033803906) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62658937) q[3];
sx q[3];
rz(-1.2703151) q[3];
sx q[3];
rz(1.8412875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5358676) q[2];
sx q[2];
rz(-2.2650227) q[2];
sx q[2];
rz(-0.65315872) q[2];
rz(-0.43290916) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(-2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2046278) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(0.2247819) q[0];
rz(-0.79554355) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(-1.068211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7601677) q[0];
sx q[0];
rz(-1.0720709) q[0];
sx q[0];
rz(2.4295761) q[0];
rz(2.3722367) q[2];
sx q[2];
rz(-2.0343668) q[2];
sx q[2];
rz(0.052841436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33419106) q[1];
sx q[1];
rz(-2.4535123) q[1];
sx q[1];
rz(-3.1013558) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8269038) q[3];
sx q[3];
rz(-2.912622) q[3];
sx q[3];
rz(0.010502149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3855359) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(-2.5059911) q[2];
rz(-0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8024837) q[0];
sx q[0];
rz(-3.1153296) q[0];
sx q[0];
rz(0.12301692) q[0];
rz(1.7440375) q[1];
sx q[1];
rz(-1.7536283) q[1];
sx q[1];
rz(-2.3618598) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5194979) q[0];
sx q[0];
rz(-2.4197398) q[0];
sx q[0];
rz(1.1546385) q[0];
rz(-pi) q[1];
rz(-0.68569195) q[2];
sx q[2];
rz(-1.4567516) q[2];
sx q[2];
rz(1.6673078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2998979) q[1];
sx q[1];
rz(-1.3720241) q[1];
sx q[1];
rz(-1.9511306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6980096) q[3];
sx q[3];
rz(-1.6774584) q[3];
sx q[3];
rz(-2.9677109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66703779) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(-3.0511268) q[2];
rz(0.41680923) q[3];
sx q[3];
rz(-2.3495245) q[3];
sx q[3];
rz(-0.99378234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823293) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(3.0521159) q[0];
rz(2.3327475) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(-2.083875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9442975) q[0];
sx q[0];
rz(-0.57387251) q[0];
sx q[0];
rz(-0.7446592) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5616778) q[2];
sx q[2];
rz(-0.55202019) q[2];
sx q[2];
rz(2.7087351) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32409975) q[1];
sx q[1];
rz(-1.9655394) q[1];
sx q[1];
rz(-0.38112074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5461033) q[3];
sx q[3];
rz(-2.3273483) q[3];
sx q[3];
rz(0.3972185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7384501) q[2];
sx q[2];
rz(-2.6046627) q[2];
sx q[2];
rz(0.37330791) q[2];
rz(0.25660723) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(3.0671425) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2557209) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(-1.4295084) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-3.1222432) q[2];
sx q[2];
rz(-0.45889284) q[2];
sx q[2];
rz(1.5677551) q[2];
rz(0.96961602) q[3];
sx q[3];
rz(-1.6027228) q[3];
sx q[3];
rz(3.0016196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
