OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(0.026260016) q[0];
rz(1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(1.4456911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64695839) q[1];
sx q[1];
rz(-2.0152433) q[1];
sx q[1];
rz(-1.7723945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20538879) q[3];
sx q[3];
rz(-0.60384149) q[3];
sx q[3];
rz(2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61629907) q[0];
sx q[0];
rz(-1.3536317) q[0];
sx q[0];
rz(2.8263682) q[0];
x q[1];
rz(1.9081406) q[2];
sx q[2];
rz(-1.0609027) q[2];
sx q[2];
rz(0.22491977) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6444781) q[1];
sx q[1];
rz(-1.9299572) q[1];
sx q[1];
rz(0.70980806) q[1];
rz(0.52821918) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(-1.5147097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(0.59282747) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.3735501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93689504) q[0];
sx q[0];
rz(-2.097313) q[0];
sx q[0];
rz(1.5419457) q[0];
x q[1];
rz(-0.85767304) q[2];
sx q[2];
rz(-1.3366941) q[2];
sx q[2];
rz(0.26299325) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2557345) q[1];
sx q[1];
rz(-2.246937) q[1];
sx q[1];
rz(2.3036792) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8009637) q[3];
sx q[3];
rz(-1.7880511) q[3];
sx q[3];
rz(1.0297071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85128879) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(1.6800866) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1876112) q[2];
sx q[2];
rz(-1.9028579) q[2];
sx q[2];
rz(-0.039475723) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(-2.2799904) q[1];
x q[2];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26966306) q[0];
sx q[0];
rz(-2.4047244) q[0];
sx q[0];
rz(2.065573) q[0];
rz(-pi) q[1];
rz(1.0079185) q[2];
sx q[2];
rz(-0.067069947) q[2];
sx q[2];
rz(-1.1941393) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.207706) q[1];
sx q[1];
rz(-1.1983786) q[1];
sx q[1];
rz(-1.9718584) q[1];
rz(0.89574121) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(0.18273396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076684549) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(-0.48318133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(-1.900308) q[0];
rz(1.5051571) q[2];
sx q[2];
rz(-2.6780431) q[2];
sx q[2];
rz(2.133873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3997765) q[1];
sx q[1];
rz(-0.92405926) q[1];
sx q[1];
rz(0.6439376) q[1];
rz(-pi) q[2];
rz(3.0275434) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(-1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-2.8894997) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(-2.8004942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3605395) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(2.1512254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7455709) q[1];
sx q[1];
rz(-0.77116291) q[1];
sx q[1];
rz(1.039617) q[1];
x q[2];
rz(0.37766405) q[3];
sx q[3];
rz(-2.0026243) q[3];
sx q[3];
rz(0.8197195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
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
rz(2.7082108) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.3716912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28384128) q[0];
sx q[0];
rz(-2.5159266) q[0];
sx q[0];
rz(-1.4197423) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3805982) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(-0.93190565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9904069) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(2.5469261) q[1];
rz(-pi) q[2];
rz(-0.38512226) q[3];
sx q[3];
rz(-2.754854) q[3];
sx q[3];
rz(0.56798565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-2.871002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71352495) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(2.9146951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2241237) q[2];
sx q[2];
rz(-1.7520095) q[2];
sx q[2];
rz(-0.62435645) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2968263) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(-1.574135) q[1];
rz(-2.050428) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(-0.8271715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.489893) q[2];
rz(2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.1606914) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(0.94296304) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-0.74238366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61481793) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(-0.74032289) q[0];
rz(0.43138357) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(-1.0647578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66354499) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(-2.7958109) q[1];
rz(-pi) q[2];
rz(-0.14245716) q[3];
sx q[3];
rz(-1.5420621) q[3];
sx q[3];
rz(1.2924259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.7998981) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-1.9146391) q[2];
sx q[2];
rz(-2.362102) q[2];
sx q[2];
rz(-1.4951928) q[2];
rz(0.7704173) q[3];
sx q[3];
rz(-0.19376783) q[3];
sx q[3];
rz(2.7289058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
