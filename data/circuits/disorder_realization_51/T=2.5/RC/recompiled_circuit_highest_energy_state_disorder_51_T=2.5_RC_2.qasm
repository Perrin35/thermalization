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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(2.0928535) q[0];
rz(0.41930786) q[1];
sx q[1];
rz(-1.3860621) q[1];
sx q[1];
rz(0.11828932) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2631419) q[0];
sx q[0];
rz(-2.1241444) q[0];
sx q[0];
rz(-0.31991495) q[0];
rz(-3.0906244) q[2];
sx q[2];
rz(-1.5870924) q[2];
sx q[2];
rz(-0.38557926) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1227545) q[1];
sx q[1];
rz(-1.600482) q[1];
sx q[1];
rz(-1.5537986) q[1];
rz(-pi) q[2];
rz(1.6476326) q[3];
sx q[3];
rz(-1.2682524) q[3];
sx q[3];
rz(1.8267814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2179541) q[2];
sx q[2];
rz(-3.1380234) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(-2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639975) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(-0.98597041) q[0];
rz(1.5892971) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(-1.5812965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88402692) q[0];
sx q[0];
rz(-2.0456893) q[0];
sx q[0];
rz(0.43254655) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5471117) q[2];
sx q[2];
rz(-1.5882601) q[2];
sx q[2];
rz(-1.6169777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1007019) q[1];
sx q[1];
rz(-1.2110079) q[1];
sx q[1];
rz(-1.3967819) q[1];
rz(-2.2038764) q[3];
sx q[3];
rz(-1.4707397) q[3];
sx q[3];
rz(1.8824739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(-1.0349234) q[2];
rz(2.6400635) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(-0.97626221) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72306776) q[0];
sx q[0];
rz(-1.1534961) q[0];
sx q[0];
rz(-1.9368517) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(-3.0012896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7492964) q[0];
sx q[0];
rz(-0.97784144) q[0];
sx q[0];
rz(-0.31813527) q[0];
rz(-pi) q[1];
rz(-0.84425462) q[2];
sx q[2];
rz(-1.2469957) q[2];
sx q[2];
rz(-2.5716788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5922484) q[1];
sx q[1];
rz(-2.1533433) q[1];
sx q[1];
rz(1.2986533) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2292333) q[3];
sx q[3];
rz(-1.7031125) q[3];
sx q[3];
rz(-1.4190471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3673765) q[2];
sx q[2];
rz(-2.1420631) q[2];
sx q[2];
rz(-2.9191169) q[2];
rz(-0.065464822) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-2.2953667) q[3];
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
rz(-1.9533933) q[0];
sx q[0];
rz(-0.024217483) q[0];
sx q[0];
rz(-0.54704332) q[0];
rz(0.25829265) q[1];
sx q[1];
rz(-0.02198418) q[1];
sx q[1];
rz(-2.8000854) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7370626) q[0];
sx q[0];
rz(-1.5642691) q[0];
sx q[0];
rz(0.0013690283) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33880587) q[2];
sx q[2];
rz(-1.1735386) q[2];
sx q[2];
rz(-0.65639979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58474135) q[1];
sx q[1];
rz(-1.9042174) q[1];
sx q[1];
rz(0.78704301) q[1];
x q[2];
rz(2.0773402) q[3];
sx q[3];
rz(-1.1501125) q[3];
sx q[3];
rz(1.3468605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7121938) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(-0.94924259) q[2];
rz(0.74887577) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6147989) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(1.5507966) q[0];
rz(1.3560449) q[1];
sx q[1];
rz(-3.1372012) q[1];
sx q[1];
rz(3.0785676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9585091) q[0];
sx q[0];
rz(-1.6367854) q[0];
sx q[0];
rz(1.6011168) q[0];
rz(2.3895415) q[2];
sx q[2];
rz(-1.1630485) q[2];
sx q[2];
rz(-0.69139451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5668727) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(1.779056) q[1];
rz(-pi) q[2];
rz(0.47565094) q[3];
sx q[3];
rz(-0.84053549) q[3];
sx q[3];
rz(-2.2851839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65681347) q[2];
sx q[2];
rz(-1.3098837) q[2];
sx q[2];
rz(-2.4978034) q[2];
rz(0.8748318) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(0.8005825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0463882) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(2.7860506) q[0];
rz(0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(0.14828646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18011151) q[0];
sx q[0];
rz(-1.5683055) q[0];
sx q[0];
rz(-2.9585696) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9695187) q[2];
sx q[2];
rz(-0.41671696) q[2];
sx q[2];
rz(-0.91815776) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10210382) q[1];
sx q[1];
rz(-1.436427) q[1];
sx q[1];
rz(0.23535442) q[1];
x q[2];
rz(2.9680524) q[3];
sx q[3];
rz(-2.1719912) q[3];
sx q[3];
rz(-1.7752826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4699012) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(3.0174603) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-2.997213) q[3];
sx q[3];
rz(0.1709443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588876) q[0];
sx q[0];
rz(-3.0165065) q[0];
sx q[0];
rz(0.74147725) q[0];
rz(-0.28400907) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(-0.31518087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717424) q[0];
sx q[0];
rz(-3.0707504) q[0];
sx q[0];
rz(1.1819345) q[0];
rz(-pi) q[1];
rz(0.20756794) q[2];
sx q[2];
rz(-1.1313442) q[2];
sx q[2];
rz(-2.07711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1629535) q[1];
sx q[1];
rz(-2.1617315) q[1];
sx q[1];
rz(1.3798452) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1180598) q[3];
sx q[3];
rz(-1.7943903) q[3];
sx q[3];
rz(1.6361039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2726941) q[2];
sx q[2];
rz(-1.0947451) q[2];
sx q[2];
rz(2.4197253) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.1451274) q[3];
sx q[3];
rz(2.0731879) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5751936) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(1.5665293) q[0];
rz(-2.9381835) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(2.4967616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5400019) q[0];
sx q[0];
rz(-0.20769697) q[0];
sx q[0];
rz(-1.5303485) q[0];
x q[1];
rz(1.5459841) q[2];
sx q[2];
rz(-0.96867079) q[2];
sx q[2];
rz(2.7393722) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1496274) q[1];
sx q[1];
rz(-2.1081165) q[1];
sx q[1];
rz(-3.0558636) q[1];
x q[2];
rz(2.6504114) q[3];
sx q[3];
rz(-1.9767663) q[3];
sx q[3];
rz(2.5570359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5975534) q[2];
sx q[2];
rz(-0.35352239) q[2];
sx q[2];
rz(2.7618347) q[2];
rz(-1.0455658) q[3];
sx q[3];
rz(-1.9087722) q[3];
sx q[3];
rz(-1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3658635) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(1.7849543) q[0];
rz(0.44048539) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(-2.4408565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2663854) q[0];
sx q[0];
rz(-0.91910494) q[0];
sx q[0];
rz(-0.65623063) q[0];
x q[1];
rz(1.3323302) q[2];
sx q[2];
rz(-0.70435134) q[2];
sx q[2];
rz(0.65504247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4266738) q[1];
sx q[1];
rz(-0.85619421) q[1];
sx q[1];
rz(-1.8780519) q[1];
rz(-2.7089617) q[3];
sx q[3];
rz(-1.5670766) q[3];
sx q[3];
rz(1.5100117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7880154) q[2];
sx q[2];
rz(-2.7699296) q[2];
sx q[2];
rz(-1.8549982) q[2];
rz(-2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(-1.9193468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-0.33682987) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754691) q[0];
sx q[0];
rz(-1.2531452) q[0];
sx q[0];
rz(-1.5730412) q[0];
x q[1];
rz(-2.1977399) q[2];
sx q[2];
rz(-2.7991931) q[2];
sx q[2];
rz(1.2417718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5975128) q[1];
sx q[1];
rz(-1.4855396) q[1];
sx q[1];
rz(-1.6823177) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1317066) q[3];
sx q[3];
rz(-1.2133382) q[3];
sx q[3];
rz(-1.4394898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6280262) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(2.8619518) q[2];
rz(0.63129342) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(-0.62425557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8414128) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(2.367876) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(2.28188) q[2];
sx q[2];
rz(-1.0621241) q[2];
sx q[2];
rz(-0.56806628) q[2];
rz(-1.5480883) q[3];
sx q[3];
rz(-1.1946214) q[3];
sx q[3];
rz(-0.0047636845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
