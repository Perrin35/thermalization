OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0041381) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(-0.11560346) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34420867) q[2];
sx q[2];
rz(-1.95103) q[2];
sx q[2];
rz(-2.3772079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73218988) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(0.7440872) q[1];
rz(-pi) q[2];
rz(2.3832541) q[3];
sx q[3];
rz(-1.7539277) q[3];
sx q[3];
rz(1.6299712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(0.83591998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0854557) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(1.4319112) q[0];
x q[1];
rz(-1.0216653) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(0.82565386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6902496) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(-2.8020225) q[1];
rz(-pi) q[2];
rz(0.42628916) q[3];
sx q[3];
rz(-2.3549035) q[3];
sx q[3];
rz(-0.28144893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(0.49697044) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.9127649) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29262221) q[0];
sx q[0];
rz(-1.6573574) q[0];
sx q[0];
rz(1.863088) q[0];
rz(-2.897981) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(-1.5872019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25698173) q[1];
sx q[1];
rz(-1.7092488) q[1];
sx q[1];
rz(0.42582663) q[1];
rz(-pi) q[2];
rz(-3.0813619) q[3];
sx q[3];
rz(-2.0647991) q[3];
sx q[3];
rz(-1.6935108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(-1.9034889) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0606196) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(0.087555126) q[0];
rz(-pi) q[1];
rz(-1.26881) q[2];
sx q[2];
rz(-0.75690818) q[2];
sx q[2];
rz(-0.41727558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31186715) q[1];
sx q[1];
rz(-2.3448181) q[1];
sx q[1];
rz(-1.7142678) q[1];
x q[2];
rz(0.55320923) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(-2.7943484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-1.0106687) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(-0.97250485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8079677) q[0];
sx q[0];
rz(-1.6728314) q[0];
sx q[0];
rz(-1.5029961) q[0];
x q[1];
rz(0.33072492) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(-1.6387788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7079066) q[1];
sx q[1];
rz(-1.3771332) q[1];
sx q[1];
rz(-2.99919) q[1];
rz(-pi) q[2];
x q[2];
rz(2.186741) q[3];
sx q[3];
rz(-1.5408437) q[3];
sx q[3];
rz(2.0314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(1.548432) q[2];
rz(1.3657773) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(-2.1877066) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71972972) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(0.37429601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4842589) q[0];
sx q[0];
rz(-1.9801635) q[0];
sx q[0];
rz(-1.9166458) q[0];
rz(-pi) q[1];
rz(1.8315973) q[2];
sx q[2];
rz(-2.0506095) q[2];
sx q[2];
rz(-0.11670437) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7172076) q[1];
sx q[1];
rz(-1.3076412) q[1];
sx q[1];
rz(-2.5724263) q[1];
rz(2.3338942) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(-0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763057) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.8315171) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6543286) q[0];
sx q[0];
rz(-1.9582821) q[0];
sx q[0];
rz(-1.5154454) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37808772) q[2];
sx q[2];
rz(-0.74947658) q[2];
sx q[2];
rz(-0.34989244) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76445626) q[1];
sx q[1];
rz(-0.38696445) q[1];
sx q[1];
rz(-2.6167234) q[1];
x q[2];
rz(0.55540107) q[3];
sx q[3];
rz(-1.5474833) q[3];
sx q[3];
rz(2.5454552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(1.6581992) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(0.21959198) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(-2.2917152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726639) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(-2.9805141) q[0];
rz(-3.0476961) q[2];
sx q[2];
rz(-2.2517423) q[2];
sx q[2];
rz(0.64881334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5596458) q[1];
sx q[1];
rz(-1.6053182) q[1];
sx q[1];
rz(2.3318961) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8306584) q[3];
sx q[3];
rz(-1.4514918) q[3];
sx q[3];
rz(2.3343149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59051096) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(-1.3809416) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.2776432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9846839) q[0];
sx q[0];
rz(-2.589993) q[0];
sx q[0];
rz(0.44996913) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4503166) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(-0.28011766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35754044) q[1];
sx q[1];
rz(-1.5338384) q[1];
sx q[1];
rz(-1.1251015) q[1];
x q[2];
rz(2.3074179) q[3];
sx q[3];
rz(-1.2779543) q[3];
sx q[3];
rz(-2.1017696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(0.35783106) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.91838592) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(2.9311438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2699379) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(-0.35738118) q[0];
x q[1];
rz(2.9774882) q[2];
sx q[2];
rz(-1.1640942) q[2];
sx q[2];
rz(1.9889529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75281843) q[1];
sx q[1];
rz(-2.7669199) q[1];
sx q[1];
rz(2.545536) q[1];
rz(-pi) q[2];
rz(-2.3841303) q[3];
sx q[3];
rz(-0.33843455) q[3];
sx q[3];
rz(-2.7528742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9615053) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5719941) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(0.16602892) q[2];
sx q[2];
rz(-2.0353073) q[2];
sx q[2];
rz(-1.4370949) q[2];
rz(-2.7995085) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
