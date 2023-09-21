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
rz(-2.3614376) q[1];
sx q[1];
rz(-1.0649788) q[1];
sx q[1];
rz(0.87632626) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1286368) q[0];
sx q[0];
rz(-2.6814333) q[0];
sx q[0];
rz(-0.53928661) q[0];
rz(0.9853739) q[2];
sx q[2];
rz(-1.4101763) q[2];
sx q[2];
rz(-2.0656245) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7180011) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(1.524339) q[1];
rz(-1.417744) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.306863) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912936) q[0];
sx q[0];
rz(-1.8739788) q[0];
sx q[0];
rz(1.1217872) q[0];
rz(-pi) q[1];
rz(2.9875056) q[2];
sx q[2];
rz(-1.0006957) q[2];
sx q[2];
rz(2.4947583) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9718711) q[1];
sx q[1];
rz(-1.1217146) q[1];
sx q[1];
rz(2.4596618) q[1];
x q[2];
rz(-1.7631093) q[3];
sx q[3];
rz(-2.8674539) q[3];
sx q[3];
rz(2.6499234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-1.0748192) q[0];
rz(0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(0.39594617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2303282) q[0];
sx q[0];
rz(-1.0697782) q[0];
sx q[0];
rz(0.45341861) q[0];
x q[1];
rz(1.2767407) q[2];
sx q[2];
rz(-1.5842423) q[2];
sx q[2];
rz(-2.3384192) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.695895) q[1];
sx q[1];
rz(-0.88765111) q[1];
sx q[1];
rz(1.5831069) q[1];
rz(-pi) q[2];
rz(1.5474165) q[3];
sx q[3];
rz(-1.3153207) q[3];
sx q[3];
rz(-1.8332831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6039156) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(2.9023857) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.061302) q[0];
sx q[0];
rz(-1.4740605) q[0];
sx q[0];
rz(-0.060782766) q[0];
x q[1];
rz(-2.4993863) q[2];
sx q[2];
rz(-1.3807447) q[2];
sx q[2];
rz(1.1914636) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4938426) q[1];
sx q[1];
rz(-2.2550681) q[1];
sx q[1];
rz(-1.1286331) q[1];
rz(-pi) q[2];
rz(-0.38846429) q[3];
sx q[3];
rz(-2.3187175) q[3];
sx q[3];
rz(-0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(-0.22339544) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-0.76104004) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2634537) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(3.1104452) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9391187) q[2];
sx q[2];
rz(-1.7624117) q[2];
sx q[2];
rz(-1.3445878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3984822) q[1];
sx q[1];
rz(-0.32144904) q[1];
sx q[1];
rz(0.076996315) q[1];
rz(-pi) q[2];
rz(1.2916336) q[3];
sx q[3];
rz(-1.2984635) q[3];
sx q[3];
rz(2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(0.87310711) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034176055) q[0];
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
rz(-2.7518227) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1205475) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(-2.8403776) q[1];
rz(-pi) q[2];
rz(1.553922) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(2.8942787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(-0.54106075) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(-3.022335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0605474) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(-2.413785) q[0];
rz(-pi) q[1];
rz(-2.1642045) q[2];
sx q[2];
rz(-2.9466629) q[2];
sx q[2];
rz(2.5755663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6357928) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(-0.15091166) q[1];
x q[2];
rz(0.70763208) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(1.5054782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83157241) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(0.88561052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53997707) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(-0.45666306) q[0];
x q[1];
rz(-1.5947072) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(-1.5469345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0091128) q[1];
sx q[1];
rz(-1.1785893) q[1];
sx q[1];
rz(-1.023804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76836821) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(2.2636569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1866859) q[0];
sx q[0];
rz(-2.1011155) q[0];
sx q[0];
rz(-0.19219877) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(-1.4093902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62577265) q[1];
sx q[1];
rz(-1.9966218) q[1];
sx q[1];
rz(1.6243402) q[1];
rz(-pi) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(2.9516454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.1970253) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4599633) q[0];
sx q[0];
rz(-2.3869793) q[0];
sx q[0];
rz(-0.90453903) q[0];
x q[1];
rz(-1.808666) q[2];
sx q[2];
rz(-1.6878205) q[2];
sx q[2];
rz(0.34840096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3665109) q[1];
sx q[1];
rz(-1.0771891) q[1];
sx q[1];
rz(-0.8702741) q[1];
rz(-pi) q[2];
rz(1.193612) q[3];
sx q[3];
rz(-1.0918655) q[3];
sx q[3];
rz(0.14648986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(0.52195436) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-0.75688731) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
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
