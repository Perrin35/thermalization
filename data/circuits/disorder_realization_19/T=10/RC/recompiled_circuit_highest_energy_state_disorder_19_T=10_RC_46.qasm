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
rz(-2.860478) q[0];
sx q[0];
rz(-1.7908362) q[0];
sx q[0];
rz(-2.1909292) q[0];
rz(0.22275337) q[1];
sx q[1];
rz(-2.8877701) q[1];
sx q[1];
rz(2.4737127) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0216103) q[0];
sx q[0];
rz(-1.9304781) q[0];
sx q[0];
rz(-1.5403454) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4761042) q[2];
sx q[2];
rz(-0.166278) q[2];
sx q[2];
rz(2.4537697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6618297) q[1];
sx q[1];
rz(-0.15108988) q[1];
sx q[1];
rz(2.2041049) q[1];
rz(-0.018785211) q[3];
sx q[3];
rz(-0.21078706) q[3];
sx q[3];
rz(-0.8715521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.604598) q[2];
sx q[2];
rz(-0.94702417) q[2];
sx q[2];
rz(0.14744082) q[2];
rz(-1.3747181) q[3];
sx q[3];
rz(-1.3876029) q[3];
sx q[3];
rz(-1.714777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0680256) q[0];
sx q[0];
rz(-1.0140714) q[0];
sx q[0];
rz(-2.9378939) q[0];
rz(-3.0286466) q[1];
sx q[1];
rz(-0.68165439) q[1];
sx q[1];
rz(3.122094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9519913) q[0];
sx q[0];
rz(-0.87965779) q[0];
sx q[0];
rz(0.36380542) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9258397) q[2];
sx q[2];
rz(-0.40242919) q[2];
sx q[2];
rz(-2.7192164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.353189) q[1];
sx q[1];
rz(-0.63523294) q[1];
sx q[1];
rz(0.53643463) q[1];
rz(-pi) q[2];
rz(0.56657378) q[3];
sx q[3];
rz(-2.1687897) q[3];
sx q[3];
rz(-1.2815042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6273177) q[2];
sx q[2];
rz(-1.5655727) q[2];
sx q[2];
rz(-2.9493098) q[2];
rz(-0.2187885) q[3];
sx q[3];
rz(-1.3008806) q[3];
sx q[3];
rz(1.7910262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.669765) q[0];
sx q[0];
rz(-0.28194031) q[0];
sx q[0];
rz(-2.6523253) q[0];
rz(-1.4504112) q[1];
sx q[1];
rz(-1.1510808) q[1];
sx q[1];
rz(-2.6077008) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.691026) q[0];
sx q[0];
rz(-1.2493142) q[0];
sx q[0];
rz(1.5430928) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0025173) q[2];
sx q[2];
rz(-1.6148881) q[2];
sx q[2];
rz(1.2033552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.008153) q[1];
sx q[1];
rz(-2.714909) q[1];
sx q[1];
rz(1.918479) q[1];
x q[2];
rz(-0.72559128) q[3];
sx q[3];
rz(-1.707336) q[3];
sx q[3];
rz(-1.8488334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78407946) q[2];
sx q[2];
rz(-2.1046905) q[2];
sx q[2];
rz(2.4927523) q[2];
rz(0.31323788) q[3];
sx q[3];
rz(-0.99010885) q[3];
sx q[3];
rz(1.8296957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.677815) q[0];
sx q[0];
rz(-0.86939335) q[0];
sx q[0];
rz(2.38548) q[0];
rz(0.4568049) q[1];
sx q[1];
rz(-2.017338) q[1];
sx q[1];
rz(2.4809428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27460262) q[0];
sx q[0];
rz(-1.1701692) q[0];
sx q[0];
rz(1.3763877) q[0];
rz(-pi) q[1];
rz(-1.2006423) q[2];
sx q[2];
rz(-1.3602019) q[2];
sx q[2];
rz(-0.38773197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7739539) q[1];
sx q[1];
rz(-0.74450297) q[1];
sx q[1];
rz(1.3013233) q[1];
rz(-pi) q[2];
rz(0.33922555) q[3];
sx q[3];
rz(-1.1407099) q[3];
sx q[3];
rz(-3.0292689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75055355) q[2];
sx q[2];
rz(-1.6497352) q[2];
sx q[2];
rz(0.98461119) q[2];
rz(-1.6330968) q[3];
sx q[3];
rz(-1.0909785) q[3];
sx q[3];
rz(2.7896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9462747) q[0];
sx q[0];
rz(-1.0349422) q[0];
sx q[0];
rz(-0.6629194) q[0];
rz(1.7341057) q[1];
sx q[1];
rz(-2.2588142) q[1];
sx q[1];
rz(-0.9443121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2812719) q[0];
sx q[0];
rz(-0.68211918) q[0];
sx q[0];
rz(-0.25231326) q[0];
rz(-2.3572982) q[2];
sx q[2];
rz(-0.81740618) q[2];
sx q[2];
rz(2.9767175) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27054003) q[1];
sx q[1];
rz(-2.4118198) q[1];
sx q[1];
rz(2.5403028) q[1];
rz(-pi) q[2];
rz(0.38864079) q[3];
sx q[3];
rz(-2.5699617) q[3];
sx q[3];
rz(2.0045668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0974836) q[2];
sx q[2];
rz(-2.5554843) q[2];
sx q[2];
rz(0.07494542) q[2];
rz(-0.99503851) q[3];
sx q[3];
rz(-1.1040265) q[3];
sx q[3];
rz(-1.1546571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.5465882) q[0];
sx q[0];
rz(-1.1052479) q[0];
sx q[0];
rz(-2.8357491) q[0];
rz(0.88286895) q[1];
sx q[1];
rz(-0.63778937) q[1];
sx q[1];
rz(0.84552228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54613333) q[0];
sx q[0];
rz(-0.89160486) q[0];
sx q[0];
rz(0.29588411) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5668198) q[2];
sx q[2];
rz(-2.8573244) q[2];
sx q[2];
rz(-2.6109744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2357422) q[1];
sx q[1];
rz(-1.0585856) q[1];
sx q[1];
rz(-0.74899574) q[1];
rz(-2.9560231) q[3];
sx q[3];
rz(-2.2757832) q[3];
sx q[3];
rz(1.5920816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4800097) q[2];
sx q[2];
rz(-2.2955387) q[2];
sx q[2];
rz(2.0704849) q[2];
rz(-0.24460159) q[3];
sx q[3];
rz(-0.94541234) q[3];
sx q[3];
rz(-1.6036114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59413183) q[0];
sx q[0];
rz(-2.5762711) q[0];
sx q[0];
rz(-2.2070337) q[0];
rz(-2.4681828) q[1];
sx q[1];
rz(-2.5118561) q[1];
sx q[1];
rz(3.0455132) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.729674) q[0];
sx q[0];
rz(-3.0239377) q[0];
sx q[0];
rz(1.3484701) q[0];
rz(2.962467) q[2];
sx q[2];
rz(-1.747588) q[2];
sx q[2];
rz(-0.13265935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9597577) q[1];
sx q[1];
rz(-1.3980306) q[1];
sx q[1];
rz(-0.29675014) q[1];
x q[2];
rz(-3.0453048) q[3];
sx q[3];
rz(-0.91595338) q[3];
sx q[3];
rz(1.8055467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7518294) q[2];
sx q[2];
rz(-1.4233754) q[2];
sx q[2];
rz(0.54330379) q[2];
rz(3.1144888) q[3];
sx q[3];
rz(-2.6684941) q[3];
sx q[3];
rz(1.0075587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2946224) q[0];
sx q[0];
rz(-0.68055081) q[0];
sx q[0];
rz(-2.3934225) q[0];
rz(-1.6698042) q[1];
sx q[1];
rz(-1.4133778) q[1];
sx q[1];
rz(2.4459623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3826806) q[0];
sx q[0];
rz(-1.6864538) q[0];
sx q[0];
rz(1.3124491) q[0];
rz(0.28098051) q[2];
sx q[2];
rz(-2.4644445) q[2];
sx q[2];
rz(-2.3317762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.66113093) q[1];
sx q[1];
rz(-0.5782402) q[1];
sx q[1];
rz(-0.64774143) q[1];
rz(-pi) q[2];
rz(-0.78277709) q[3];
sx q[3];
rz(-0.15412384) q[3];
sx q[3];
rz(-1.085404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2601605) q[2];
sx q[2];
rz(-2.3023534) q[2];
sx q[2];
rz(-0.053675573) q[2];
rz(1.7565049) q[3];
sx q[3];
rz(-1.3421007) q[3];
sx q[3];
rz(-2.9746941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3593035) q[0];
sx q[0];
rz(-1.8817236) q[0];
sx q[0];
rz(-0.89349973) q[0];
rz(-2.9529849) q[1];
sx q[1];
rz(-0.74917561) q[1];
sx q[1];
rz(-0.20208727) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2938439) q[0];
sx q[0];
rz(-1.6007152) q[0];
sx q[0];
rz(-0.82383967) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0918504) q[2];
sx q[2];
rz(-2.0936077) q[2];
sx q[2];
rz(-0.78429121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1043865) q[1];
sx q[1];
rz(-1.361711) q[1];
sx q[1];
rz(-2.1165743) q[1];
rz(-pi) q[2];
rz(-1.8340183) q[3];
sx q[3];
rz(-1.211245) q[3];
sx q[3];
rz(1.5375801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9933652) q[2];
sx q[2];
rz(-1.8348285) q[2];
sx q[2];
rz(0.6692878) q[2];
rz(0.78684849) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(0.70080668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65449077) q[0];
sx q[0];
rz(-0.50906068) q[0];
sx q[0];
rz(1.2855726) q[0];
rz(2.9945943) q[1];
sx q[1];
rz(-1.1451984) q[1];
sx q[1];
rz(2.1910892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5326544) q[0];
sx q[0];
rz(-1.0883347) q[0];
sx q[0];
rz(2.5368693) q[0];
rz(1.6754402) q[2];
sx q[2];
rz(-0.54000914) q[2];
sx q[2];
rz(-1.4252942) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1482137) q[1];
sx q[1];
rz(-1.0328173) q[1];
sx q[1];
rz(-0.64126539) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13895039) q[3];
sx q[3];
rz(-1.6873056) q[3];
sx q[3];
rz(-2.9187035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6805083) q[2];
sx q[2];
rz(-0.66404873) q[2];
sx q[2];
rz(2.7456679) q[2];
rz(2.8094273) q[3];
sx q[3];
rz(-1.1484523) q[3];
sx q[3];
rz(-0.90126669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.111515) q[0];
sx q[0];
rz(-0.60965309) q[0];
sx q[0];
rz(-1.6994221) q[0];
rz(0.036399966) q[1];
sx q[1];
rz(-1.8489238) q[1];
sx q[1];
rz(1.363149) q[1];
rz(-1.8698975) q[2];
sx q[2];
rz(-2.3334685) q[2];
sx q[2];
rz(0.12167385) q[2];
rz(2.4717109) q[3];
sx q[3];
rz(-1.6139779) q[3];
sx q[3];
rz(-1.1077751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
