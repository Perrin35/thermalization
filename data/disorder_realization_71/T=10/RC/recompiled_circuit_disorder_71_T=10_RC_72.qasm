OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(2.2652664) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0759461) q[0];
sx q[0];
rz(-1.8008721) q[0];
sx q[0];
rz(0.40212698) q[0];
rz(-0.19199065) q[2];
sx q[2];
rz(-0.99388323) q[2];
sx q[2];
rz(-2.752395) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0383366) q[1];
sx q[1];
rz(-1.5557319) q[1];
sx q[1];
rz(-1.9008093) q[1];
x q[2];
rz(-1.7238486) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-2.412964) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6073608) q[0];
sx q[0];
rz(-2.605654) q[0];
sx q[0];
rz(2.1952654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1463257) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(-1.0075943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.892131) q[1];
sx q[1];
rz(-0.79626894) q[1];
sx q[1];
rz(2.4888121) q[1];
rz(-pi) q[2];
rz(1.3014684) q[3];
sx q[3];
rz(-1.519031) q[3];
sx q[3];
rz(2.2477637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-0.43593105) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(2.7456465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1186737) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(2.2454717) q[0];
rz(-pi) q[1];
rz(3.1275438) q[2];
sx q[2];
rz(-1.8648246) q[2];
sx q[2];
rz(-2.369898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6763941) q[1];
sx q[1];
rz(-0.68323831) q[1];
sx q[1];
rz(0.015124358) q[1];
x q[2];
rz(1.5474165) q[3];
sx q[3];
rz(-1.8262719) q[3];
sx q[3];
rz(-1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(2.9023857) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6452091) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(1.4738826) q[0];
rz(-pi) q[1];
rz(-1.3350305) q[2];
sx q[2];
rz(-0.94199099) q[2];
sx q[2];
rz(-2.6218888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64775002) q[1];
sx q[1];
rz(-2.2550681) q[1];
sx q[1];
rz(-1.1286331) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1831746) q[3];
sx q[3];
rz(-0.82510199) q[3];
sx q[3];
rz(-2.2771698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(-0.24211611) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(0.76104004) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431664) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(-1.4545928) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9391187) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(-1.7970049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(-0.076996315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3628999) q[3];
sx q[3];
rz(-0.38749309) q[3];
sx q[3];
rz(3.1117698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(0.87310711) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(-2.81566) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1074166) q[0];
sx q[0];
rz(-2.0909485) q[0];
sx q[0];
rz(-1.1362032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44547703) q[2];
sx q[2];
rz(-1.8468841) q[2];
sx q[2];
rz(0.38976994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1833916) q[1];
sx q[1];
rz(-0.32918731) q[1];
sx q[1];
rz(0.42894657) q[1];
x q[2];
rz(-1.553922) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(-2.8942787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(2.4874172) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(0.54106075) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-0.72189271) q[0];
rz(1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(3.022335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0605474) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(0.72780769) q[0];
rz(-1.7330405) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(-0.42019368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1489361) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(2.0723144) q[1];
rz(-pi) q[2];
rz(-1.2001541) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3100202) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(0.88561052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793593) q[0];
sx q[0];
rz(-2.057312) q[0];
sx q[0];
rz(1.3079206) q[0];
x q[1];
rz(-1.5947072) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(1.5946582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9319716) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(-0.45100905) q[1];
x q[2];
rz(-0.86730154) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(-2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(2.0057604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95490676) q[0];
sx q[0];
rz(-2.1011155) q[0];
sx q[0];
rz(-2.9493939) q[0];
x q[1];
rz(-0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(-1.6795295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.218704) q[1];
sx q[1];
rz(-1.6195546) q[1];
sx q[1];
rz(2.7152275) q[1];
x q[2];
rz(-1.0335835) q[3];
sx q[3];
rz(-2.4164003) q[3];
sx q[3];
rz(2.4560526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.2783485) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.1970253) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9997864) q[0];
sx q[0];
rz(-2.1394661) q[0];
sx q[0];
rz(-2.6151711) q[0];
rz(3.0212101) q[2];
sx q[2];
rz(-1.8070081) q[2];
sx q[2];
rz(1.2506968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7750818) q[1];
sx q[1];
rz(-2.0644036) q[1];
sx q[1];
rz(0.8702741) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61693807) q[3];
sx q[3];
rz(-0.60041282) q[3];
sx q[3];
rz(-0.85655772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(0.79997921) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(2.7675046) q[2];
sx q[2];
rz(-1.2384407) q[2];
sx q[2];
rz(0.25638914) q[2];
rz(0.92267358) q[3];
sx q[3];
rz(-2.1019983) q[3];
sx q[3];
rz(-2.0207873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
