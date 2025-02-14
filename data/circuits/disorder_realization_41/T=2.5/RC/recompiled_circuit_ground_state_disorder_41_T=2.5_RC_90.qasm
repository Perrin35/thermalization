OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(-1.2394152) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83441855) q[0];
sx q[0];
rz(-2.7980013) q[0];
sx q[0];
rz(-0.8309721) q[0];
rz(2.5074717) q[2];
sx q[2];
rz(-2.8180052) q[2];
sx q[2];
rz(-0.35558082) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2762294) q[1];
sx q[1];
rz(-1.6003834) q[1];
sx q[1];
rz(0.33025708) q[1];
rz(-2.4437583) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(-0.35240155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4702845) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(-0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(-2.9644137) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.141356) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(-0.29362383) q[0];
rz(-2.1108421) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(-2.7986599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200877) q[0];
sx q[0];
rz(-1.3653127) q[0];
sx q[0];
rz(1.0387102) q[0];
x q[1];
rz(0.12983506) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(2.2463727) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56655069) q[1];
sx q[1];
rz(-0.88156619) q[1];
sx q[1];
rz(-1.8495967) q[1];
rz(-pi) q[2];
rz(-2.2003056) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(-2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0127516) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-3.1208842) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(-0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6193806) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(-0.0044599175) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(2.3562145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1566832) q[0];
sx q[0];
rz(-1.4117318) q[0];
sx q[0];
rz(1.6522264) q[0];
rz(1.8846604) q[2];
sx q[2];
rz(-1.403927) q[2];
sx q[2];
rz(0.97970574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1238034) q[1];
sx q[1];
rz(-2.5964156) q[1];
sx q[1];
rz(1.4383609) q[1];
x q[2];
rz(-0.80180577) q[3];
sx q[3];
rz(-2.2983645) q[3];
sx q[3];
rz(1.3096732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(-2.2564364) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(-2.0451827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.261026) q[0];
sx q[0];
rz(-1.9184904) q[0];
sx q[0];
rz(0.33190042) q[0];
x q[1];
rz(3.1109516) q[2];
sx q[2];
rz(-1.6302135) q[2];
sx q[2];
rz(-2.5031646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1433574) q[1];
sx q[1];
rz(-1.1117742) q[1];
sx q[1];
rz(1.8602536) q[1];
rz(-2.0799594) q[3];
sx q[3];
rz(-0.90583767) q[3];
sx q[3];
rz(3.0388415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61802822) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(-0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(0.88328254) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4516975) q[0];
sx q[0];
rz(-0.30245936) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(0.11058552) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(1.0708403) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12884451) q[0];
sx q[0];
rz(-1.8086047) q[0];
sx q[0];
rz(-1.1046011) q[0];
rz(-0.51009615) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(2.9053743) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22510281) q[1];
sx q[1];
rz(-1.9472633) q[1];
sx q[1];
rz(0.55209362) q[1];
rz(-pi) q[2];
rz(2.4399906) q[3];
sx q[3];
rz(-2.1013341) q[3];
sx q[3];
rz(-0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.11833) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.2888541) q[0];
rz(-1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.9992794) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5286251) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(0.27530833) q[0];
rz(-pi) q[1];
rz(0.062131957) q[2];
sx q[2];
rz(-1.7781864) q[2];
sx q[2];
rz(1.0203938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21396046) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(0.63251782) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8303538) q[3];
sx q[3];
rz(-1.7390307) q[3];
sx q[3];
rz(0.5420891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3219354) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(-0.22404484) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(0.45101756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7190349) q[0];
sx q[0];
rz(-1.8717878) q[0];
sx q[0];
rz(-1.6331571) q[0];
rz(1.9372378) q[2];
sx q[2];
rz(-1.7882573) q[2];
sx q[2];
rz(1.6161473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1731733) q[1];
sx q[1];
rz(-2.0035158) q[1];
sx q[1];
rz(1.7041901) q[1];
rz(0.0093098442) q[3];
sx q[3];
rz(-2.3311989) q[3];
sx q[3];
rz(1.9805465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1958127) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(1.3607402) q[2];
rz(-1.0055379) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(2.4355167) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(0.85618883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5917011) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.9338305) q[0];
rz(-pi) q[1];
rz(-2.0308308) q[2];
sx q[2];
rz(-2.3383814) q[2];
sx q[2];
rz(-0.75424657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44196196) q[1];
sx q[1];
rz(-0.61304997) q[1];
sx q[1];
rz(0.35787257) q[1];
rz(-pi) q[2];
rz(1.8185577) q[3];
sx q[3];
rz(-2.6399586) q[3];
sx q[3];
rz(-0.35546965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(-0.16743463) q[2];
rz(-1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3779959) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(-0.30002123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094291501) q[0];
sx q[0];
rz(-2.8326748) q[0];
sx q[0];
rz(0.39157097) q[0];
rz(-pi) q[1];
rz(0.062252684) q[2];
sx q[2];
rz(-2.5635898) q[2];
sx q[2];
rz(-2.5684772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2310861) q[1];
sx q[1];
rz(-2.4795737) q[1];
sx q[1];
rz(-0.20259095) q[1];
rz(3.0172273) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(1.82064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(0.023177711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4172702) q[0];
sx q[0];
rz(-0.5529772) q[0];
sx q[0];
rz(2.824845) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2397348) q[2];
sx q[2];
rz(-1.0791856) q[2];
sx q[2];
rz(2.7265446) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0129628) q[1];
sx q[1];
rz(-2.9010765) q[1];
sx q[1];
rz(0.81324767) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7816824) q[3];
sx q[3];
rz(-2.4004705) q[3];
sx q[3];
rz(1.9884584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(0.49087697) q[2];
rz(-2.0513746) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(0.31392613) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8563817) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(-2.8648227) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(2.4527373) q[2];
sx q[2];
rz(-1.3313455) q[2];
sx q[2];
rz(-0.81305885) q[2];
rz(1.5506977) q[3];
sx q[3];
rz(-2.1051959) q[3];
sx q[3];
rz(0.013230562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
