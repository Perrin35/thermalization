OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(0.63011516) q[0];
rz(-2.6684127) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(2.6602886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66274553) q[1];
sx q[1];
rz(-1.6531684) q[1];
sx q[1];
rz(-1.808951) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8880635) q[3];
sx q[3];
rz(-2.5228365) q[3];
sx q[3];
rz(1.7537102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(1.1260024) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(0.69357187) q[0];
rz(-1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-2.9512761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961401) q[0];
sx q[0];
rz(-1.7691233) q[0];
sx q[0];
rz(1.1597) q[0];
rz(-pi) q[1];
rz(-0.53493494) q[2];
sx q[2];
rz(-1.9133647) q[2];
sx q[2];
rz(-1.5869706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7029876) q[1];
sx q[1];
rz(-2.4191796) q[1];
sx q[1];
rz(1.5370675) q[1];
x q[2];
rz(-0.5204366) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4753715) q[0];
sx q[0];
rz(-0.9733805) q[0];
sx q[0];
rz(-1.4587547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9343611) q[2];
sx q[2];
rz(-1.6336294) q[2];
sx q[2];
rz(2.732423) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4856845) q[1];
sx q[1];
rz(-2.471172) q[1];
sx q[1];
rz(0.1078492) q[1];
rz(-pi) q[2];
rz(-2.2218024) q[3];
sx q[3];
rz(-1.8120822) q[3];
sx q[3];
rz(-1.9581219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766749) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(-1.2491559) q[0];
rz(-pi) q[1];
rz(-1.2816216) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(0.5493872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0399196) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(2.4458829) q[1];
rz(-pi) q[2];
rz(-1.6468871) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2356448) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(0.91919351) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910687) q[0];
sx q[0];
rz(-2.2254125) q[0];
sx q[0];
rz(3.0380681) q[0];
rz(-0.23767383) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(1.8397457) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6225699) q[1];
sx q[1];
rz(-1.4624603) q[1];
sx q[1];
rz(0.83801724) q[1];
rz(-pi) q[2];
rz(-0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(-0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-1.3549266) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(2.5700263) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.556658) q[0];
sx q[0];
rz(-2.1437763) q[0];
sx q[0];
rz(-2.1893188) q[0];
rz(-0.90494855) q[2];
sx q[2];
rz(-1.6709423) q[2];
sx q[2];
rz(0.03749321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.188365) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(0.26785775) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70116331) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(0.13247709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-0.71227658) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(0.66551048) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367046) q[0];
sx q[0];
rz(-2.8993653) q[0];
sx q[0];
rz(-1.5664943) q[0];
rz(-0.14851103) q[2];
sx q[2];
rz(-2.2627137) q[2];
sx q[2];
rz(1.4505381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.682714) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-2.0290124) q[1];
x q[2];
rz(1.3286367) q[3];
sx q[3];
rz(-1.7112964) q[3];
sx q[3];
rz(1.5797918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20913798) q[0];
sx q[0];
rz(-1.4203686) q[0];
sx q[0];
rz(-1.9183137) q[0];
x q[1];
rz(-0.4523925) q[2];
sx q[2];
rz(-2.6466742) q[2];
sx q[2];
rz(-0.53168833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.056811) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(1.0952428) q[1];
x q[2];
rz(2.268928) q[3];
sx q[3];
rz(-0.41655311) q[3];
sx q[3];
rz(-1.0636898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(2.4832446) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(0.13959612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9717279) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(1.6318897) q[0];
rz(0.058768674) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(-1.2440484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3176206) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(0.13336639) q[1];
x q[2];
rz(-0.052501909) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(-2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6167986) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-0.087879114) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.6569482) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(1.0111965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-3.0849506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6757322) q[1];
sx q[1];
rz(-1.4714186) q[1];
sx q[1];
rz(-1.9468716) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0049403355) q[3];
sx q[3];
rz(-2.2757109) q[3];
sx q[3];
rz(-2.4019394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(1.4245695) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(0.88541661) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
