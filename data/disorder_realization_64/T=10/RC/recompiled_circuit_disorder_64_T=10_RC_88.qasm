OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(3.1199772) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3992213) q[0];
sx q[0];
rz(-1.0065777) q[0];
sx q[0];
rz(-0.39682927) q[0];
rz(-pi) q[1];
rz(3.0388019) q[2];
sx q[2];
rz(-1.3445026) q[2];
sx q[2];
rz(-3.0068827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0845619) q[1];
sx q[1];
rz(-1.44094) q[1];
sx q[1];
rz(2.3741541) q[1];
x q[2];
rz(2.0256151) q[3];
sx q[3];
rz(-0.67512073) q[3];
sx q[3];
rz(-1.7636553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(2.0143051) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(-0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-3.1399472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5434108) q[0];
sx q[0];
rz(-3.0953005) q[0];
sx q[0];
rz(1.339519) q[0];
rz(-pi) q[1];
rz(1.8284945) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(-1.4305654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.590608) q[1];
sx q[1];
rz(-1.814517) q[1];
sx q[1];
rz(-2.3477712) q[1];
rz(2.1529589) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(-2.8305588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(2.0522096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82747805) q[0];
sx q[0];
rz(-2.4347322) q[0];
sx q[0];
rz(0.88721888) q[0];
rz(-pi) q[1];
rz(1.9956279) q[2];
sx q[2];
rz(-2.1263188) q[2];
sx q[2];
rz(2.5512763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5730126) q[1];
sx q[1];
rz(-1.9553361) q[1];
sx q[1];
rz(-0.09210715) q[1];
rz(-pi) q[2];
rz(-0.55604071) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(-1.3211105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(-2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.115288) q[0];
sx q[0];
rz(-1.3867154) q[0];
sx q[0];
rz(0.91362761) q[0];
rz(-pi) q[1];
rz(-0.6511351) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(-0.091094253) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4543912) q[1];
sx q[1];
rz(-1.9481716) q[1];
sx q[1];
rz(-0.45339938) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30712819) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(-3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(0.75677538) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.416052) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(2.1062772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0400378) q[0];
sx q[0];
rz(-2.3674175) q[0];
sx q[0];
rz(-0.9206307) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0683124) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(-2.5626593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6742179) q[1];
sx q[1];
rz(-1.6661223) q[1];
sx q[1];
rz(1.5159038) q[1];
rz(-0.27541311) q[3];
sx q[3];
rz(-1.3684891) q[3];
sx q[3];
rz(-1.3261258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(0.76888293) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(-1.6736354) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(0.2072269) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6672872) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(-2.5179203) q[0];
rz(1.1057165) q[2];
sx q[2];
rz(-1.0208703) q[2];
sx q[2];
rz(-1.20649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8689649) q[1];
sx q[1];
rz(-2.7902863) q[1];
sx q[1];
rz(0.79877324) q[1];
x q[2];
rz(1.7436142) q[3];
sx q[3];
rz(-0.76068766) q[3];
sx q[3];
rz(2.2067604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-2.4198789) q[2];
rz(1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(-0.2917372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7043982) q[0];
sx q[0];
rz(-1.6219553) q[0];
sx q[0];
rz(-0.41914661) q[0];
rz(-0.10642274) q[2];
sx q[2];
rz(-1.8723882) q[2];
sx q[2];
rz(-2.0472722) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8768423) q[1];
sx q[1];
rz(-1.8732338) q[1];
sx q[1];
rz(1.7251863) q[1];
rz(-1.8802059) q[3];
sx q[3];
rz(-1.0687807) q[3];
sx q[3];
rz(1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(0.062019197) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-2.4777381) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52866919) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(2.5445166) q[0];
x q[1];
rz(-0.27997048) q[2];
sx q[2];
rz(-1.0500056) q[2];
sx q[2];
rz(-2.1998646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16756032) q[1];
sx q[1];
rz(-1.4640199) q[1];
sx q[1];
rz(2.0978931) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.557425) q[3];
sx q[3];
rz(-1.468588) q[3];
sx q[3];
rz(2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(0.46554309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532928) q[0];
sx q[0];
rz(-2.4294937) q[0];
sx q[0];
rz(-3.0112991) q[0];
x q[1];
rz(-3.0131857) q[2];
sx q[2];
rz(-1.9799211) q[2];
sx q[2];
rz(2.1985334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1411966) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(2.6617906) q[1];
rz(-pi) q[2];
rz(-0.25119987) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(2.420345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.643606) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(0.75138584) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-2.5591992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2975785) q[0];
sx q[0];
rz(-1.7739002) q[0];
sx q[0];
rz(2.2556979) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2920612) q[2];
sx q[2];
rz(-2.6346452) q[2];
sx q[2];
rz(-0.35001937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(0.99517676) q[1];
rz(-pi) q[2];
rz(-0.21721812) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(0.56009968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-1.273524) q[2];
sx q[2];
rz(-2.3625629) q[2];
sx q[2];
rz(0.071803781) q[2];
rz(1.0675666) q[3];
sx q[3];
rz(-1.1451086) q[3];
sx q[3];
rz(1.2721636) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];