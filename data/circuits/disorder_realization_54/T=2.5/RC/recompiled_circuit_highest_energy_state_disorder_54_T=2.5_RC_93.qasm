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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(-0.56587926) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(4.9257435) q[1];
sx q[1];
rz(9.9590376) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256828) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(1.9749667) q[0];
rz(-2.7083394) q[2];
sx q[2];
rz(-2.5009071) q[2];
sx q[2];
rz(-2.6034466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6001662) q[1];
sx q[1];
rz(-1.2703589) q[1];
sx q[1];
rz(1.5460127) q[1];
rz(-1.322106) q[3];
sx q[3];
rz(-1.1258896) q[3];
sx q[3];
rz(1.9865195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1285706) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(0.070076076) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(-2.8299423) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79651901) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(-2.2497798) q[0];
rz(-3.0333854) q[2];
sx q[2];
rz(-1.3851067) q[2];
sx q[2];
rz(-0.78150392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5923002) q[1];
sx q[1];
rz(-0.29634297) q[1];
sx q[1];
rz(-0.94409385) q[1];
x q[2];
rz(-2.8007224) q[3];
sx q[3];
rz(-2.1580527) q[3];
sx q[3];
rz(1.0140918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0011757294) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(2.1902093) q[2];
rz(2.9133255) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(-1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(0.1618298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069050633) q[0];
sx q[0];
rz(-0.74063535) q[0];
sx q[0];
rz(0.97461755) q[0];
x q[1];
rz(2.929737) q[2];
sx q[2];
rz(-0.39688928) q[2];
sx q[2];
rz(-2.3174163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.075167478) q[1];
sx q[1];
rz(-2.4801765) q[1];
sx q[1];
rz(-0.39638806) q[1];
x q[2];
rz(-0.39997081) q[3];
sx q[3];
rz(-2.2999527) q[3];
sx q[3];
rz(3.0384881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(-0.67489135) q[0];
rz(-0.15448013) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(0.45796576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77112326) q[0];
sx q[0];
rz(-1.1827977) q[0];
sx q[0];
rz(2.6113135) q[0];
rz(-2.8343924) q[2];
sx q[2];
rz(-1.9824902) q[2];
sx q[2];
rz(1.1634942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.612466) q[1];
sx q[1];
rz(-1.3867897) q[1];
sx q[1];
rz(0.37671582) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9513957) q[3];
sx q[3];
rz(-0.51589291) q[3];
sx q[3];
rz(-0.35515337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0080879) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(2.8221455) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(0.5602079) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(-1.9200448) q[0];
rz(-0.23280652) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(-0.98308841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670416) q[0];
sx q[0];
rz(-0.75387757) q[0];
sx q[0];
rz(2.3874832) q[0];
rz(-2.1914761) q[2];
sx q[2];
rz(-1.663637) q[2];
sx q[2];
rz(-2.8454121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32199643) q[1];
sx q[1];
rz(-0.23769544) q[1];
sx q[1];
rz(2.4758215) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9748989) q[3];
sx q[3];
rz(-1.1536666) q[3];
sx q[3];
rz(-0.47458173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(1.075607) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(1.7042473) q[0];
rz(-3.0628693) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054436) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(2.8092119) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70912045) q[2];
sx q[2];
rz(-2.4599584) q[2];
sx q[2];
rz(2.1644985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50997616) q[1];
sx q[1];
rz(-2.2215034) q[1];
sx q[1];
rz(-1.0751616) q[1];
x q[2];
rz(-0.6491361) q[3];
sx q[3];
rz(-2.2723291) q[3];
sx q[3];
rz(2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(-0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5439344) q[0];
sx q[0];
rz(-1.1827396) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(-1.687382) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(1.1167663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46653173) q[0];
sx q[0];
rz(-1.2197615) q[0];
sx q[0];
rz(0.37295966) q[0];
rz(0.26979763) q[2];
sx q[2];
rz(-1.0493663) q[2];
sx q[2];
rz(0.21873378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36193725) q[1];
sx q[1];
rz(-1.9101659) q[1];
sx q[1];
rz(1.9091868) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0159166) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(-2.3762459) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(0.76559693) q[0];
rz(-0.87144026) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(-1.8188459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8894371) q[0];
sx q[0];
rz(-1.8312635) q[0];
sx q[0];
rz(-1.5171264) q[0];
rz(-pi) q[1];
rz(2.2501588) q[2];
sx q[2];
rz(-1.3855033) q[2];
sx q[2];
rz(-0.44033716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.821956) q[1];
sx q[1];
rz(-1.0624891) q[1];
sx q[1];
rz(1.5519358) q[1];
rz(-0.83006184) q[3];
sx q[3];
rz(-0.30127159) q[3];
sx q[3];
rz(-0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65779296) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(-2.9198809) q[2];
rz(1.0929557) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-0.67813412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073864989) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-0.27716032) q[0];
rz(-2.3932338) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(-1.998273) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3353334) q[0];
sx q[0];
rz(-1.8758095) q[0];
sx q[0];
rz(2.9312031) q[0];
rz(-pi) q[1];
rz(-2.2758946) q[2];
sx q[2];
rz(-1.7324466) q[2];
sx q[2];
rz(-0.18582144) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1705679) q[1];
sx q[1];
rz(-1.2538365) q[1];
sx q[1];
rz(-0.16522878) q[1];
rz(-pi) q[2];
rz(1.7420898) q[3];
sx q[3];
rz(-1.2138324) q[3];
sx q[3];
rz(-1.5689284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(-0.048967036) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5628409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9492968) q[0];
sx q[0];
rz(-1.6598633) q[0];
sx q[0];
rz(1.7013447) q[0];
x q[1];
rz(-1.3124002) q[2];
sx q[2];
rz(-2.5157197) q[2];
sx q[2];
rz(1.3385119) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.409118) q[1];
sx q[1];
rz(-0.98608398) q[1];
sx q[1];
rz(0.52701108) q[1];
rz(-pi) q[2];
rz(2.3046523) q[3];
sx q[3];
rz(-1.6815974) q[3];
sx q[3];
rz(-1.5821804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(-2.9760402) q[2];
rz(0.60756573) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7623357) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
rz(1.634585) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(2.3661939) q[2];
sx q[2];
rz(-1.4682894) q[2];
sx q[2];
rz(-0.52649047) q[2];
rz(0.60579469) q[3];
sx q[3];
rz(-1.7346564) q[3];
sx q[3];
rz(0.32614542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
