OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65741444) q[0];
sx q[0];
rz(1.2793469) q[0];
sx q[0];
rz(10.068324) q[0];
rz(-2.9930288) q[1];
sx q[1];
rz(-0.64639503) q[1];
sx q[1];
rz(0.57599154) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9849562) q[0];
sx q[0];
rz(-0.63581563) q[0];
sx q[0];
rz(-2.1106857) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0895626) q[2];
sx q[2];
rz(-0.33809755) q[2];
sx q[2];
rz(-1.372499) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.483584) q[1];
sx q[1];
rz(-1.600126) q[1];
sx q[1];
rz(1.6701677) q[1];
rz(-pi) q[2];
rz(2.9408119) q[3];
sx q[3];
rz(-2.3030278) q[3];
sx q[3];
rz(-1.4522566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6286455) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(2.4224572) q[2];
rz(0.61270815) q[3];
sx q[3];
rz(-2.6477974) q[3];
sx q[3];
rz(1.6814211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7270108) q[0];
sx q[0];
rz(-0.56986037) q[0];
sx q[0];
rz(-1.7852831) q[0];
rz(0.59056979) q[1];
sx q[1];
rz(-1.5868203) q[1];
sx q[1];
rz(-2.2812567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7588494) q[0];
sx q[0];
rz(-0.68419391) q[0];
sx q[0];
rz(-1.8608776) q[0];
rz(2.8639272) q[2];
sx q[2];
rz(-0.57304731) q[2];
sx q[2];
rz(-2.3615357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0257657) q[1];
sx q[1];
rz(-1.2526667) q[1];
sx q[1];
rz(2.8092188) q[1];
rz(-pi) q[2];
rz(-1.6708871) q[3];
sx q[3];
rz(-2.2948569) q[3];
sx q[3];
rz(1.0919746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4071953) q[2];
sx q[2];
rz(-2.037214) q[2];
sx q[2];
rz(2.8779136) q[2];
rz(-3.0401491) q[3];
sx q[3];
rz(-1.0976617) q[3];
sx q[3];
rz(1.6775848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9145255) q[0];
sx q[0];
rz(-2.3093746) q[0];
sx q[0];
rz(-1.2365923) q[0];
rz(2.6446758) q[1];
sx q[1];
rz(-2.2221815) q[1];
sx q[1];
rz(-0.20981728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1001756) q[0];
sx q[0];
rz(-0.13209535) q[0];
sx q[0];
rz(2.9939674) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0044697) q[2];
sx q[2];
rz(-1.1242766) q[2];
sx q[2];
rz(1.6046804) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1637437) q[1];
sx q[1];
rz(-1.5340549) q[1];
sx q[1];
rz(-1.0669462) q[1];
rz(2.9852226) q[3];
sx q[3];
rz(-0.84330446) q[3];
sx q[3];
rz(-2.7263977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0158374) q[2];
sx q[2];
rz(-1.1344942) q[2];
sx q[2];
rz(-2.8718359) q[2];
rz(2.6848327) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(0.30341283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.03511) q[0];
sx q[0];
rz(-0.39316097) q[0];
sx q[0];
rz(3.0495354) q[0];
rz(2.5606142) q[1];
sx q[1];
rz(-2.5216504) q[1];
sx q[1];
rz(0.43749014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96707805) q[0];
sx q[0];
rz(-2.1319492) q[0];
sx q[0];
rz(-0.18995096) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3542489) q[2];
sx q[2];
rz(-1.4581973) q[2];
sx q[2];
rz(-1.4527904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0713378) q[1];
sx q[1];
rz(-1.1769088) q[1];
sx q[1];
rz(0.27166922) q[1];
x q[2];
rz(-2.8169223) q[3];
sx q[3];
rz(-0.89327795) q[3];
sx q[3];
rz(0.087740104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1979394) q[2];
sx q[2];
rz(-1.3760171) q[2];
sx q[2];
rz(-0.1680689) q[2];
rz(-2.4456612) q[3];
sx q[3];
rz(-1.1974502) q[3];
sx q[3];
rz(-1.7456938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.593852) q[0];
sx q[0];
rz(-2.934179) q[0];
sx q[0];
rz(-0.70163027) q[0];
rz(-0.54948366) q[1];
sx q[1];
rz(-1.0797078) q[1];
sx q[1];
rz(2.2061677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.219329) q[0];
sx q[0];
rz(-0.50271243) q[0];
sx q[0];
rz(-1.034584) q[0];
rz(-pi) q[1];
rz(-1.6239426) q[2];
sx q[2];
rz(-1.4477483) q[2];
sx q[2];
rz(1.4966663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0282057) q[1];
sx q[1];
rz(-1.0019687) q[1];
sx q[1];
rz(-0.48888205) q[1];
rz(-1.6356556) q[3];
sx q[3];
rz(-1.4993414) q[3];
sx q[3];
rz(-0.91141975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51297411) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(-2.4690907) q[2];
rz(-2.3682112) q[3];
sx q[3];
rz(-1.0459817) q[3];
sx q[3];
rz(-0.33236233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32909876) q[0];
sx q[0];
rz(-2.2342873) q[0];
sx q[0];
rz(-1.4688274) q[0];
rz(2.2724197) q[1];
sx q[1];
rz(-2.4425127) q[1];
sx q[1];
rz(-0.3259784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7237968) q[0];
sx q[0];
rz(-1.4362596) q[0];
sx q[0];
rz(2.8409182) q[0];
rz(2.432304) q[2];
sx q[2];
rz(-0.69586674) q[2];
sx q[2];
rz(0.33791956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0404929) q[1];
sx q[1];
rz(-1.0940576) q[1];
sx q[1];
rz(0.65495305) q[1];
rz(-pi) q[2];
x q[2];
rz(0.066727344) q[3];
sx q[3];
rz(-1.9983833) q[3];
sx q[3];
rz(-0.59837435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.088923067) q[2];
sx q[2];
rz(-1.5648877) q[2];
sx q[2];
rz(-0.79724533) q[2];
rz(0.097660216) q[3];
sx q[3];
rz(-2.3982513) q[3];
sx q[3];
rz(1.2231479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.866211) q[0];
sx q[0];
rz(-0.99358639) q[0];
sx q[0];
rz(2.1589101) q[0];
rz(-1.6536225) q[1];
sx q[1];
rz(-2.5762312) q[1];
sx q[1];
rz(0.32378325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2764012) q[0];
sx q[0];
rz(-2.8439024) q[0];
sx q[0];
rz(-0.28930347) q[0];
rz(-1.9672438) q[2];
sx q[2];
rz(-0.65047036) q[2];
sx q[2];
rz(-2.5983368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6205755) q[1];
sx q[1];
rz(-1.4572826) q[1];
sx q[1];
rz(-0.6599636) q[1];
rz(-pi) q[2];
rz(-2.5681132) q[3];
sx q[3];
rz(-0.57807589) q[3];
sx q[3];
rz(-2.7356407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4299778) q[2];
sx q[2];
rz(-2.6275676) q[2];
sx q[2];
rz(-0.035710486) q[2];
rz(-2.4014373) q[3];
sx q[3];
rz(-1.454708) q[3];
sx q[3];
rz(1.9945701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78298727) q[0];
sx q[0];
rz(-2.7645223) q[0];
sx q[0];
rz(-2.9407035) q[0];
rz(2.5783077) q[1];
sx q[1];
rz(-1.3689684) q[1];
sx q[1];
rz(2.7422781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4056643) q[0];
sx q[0];
rz(-2.9333994) q[0];
sx q[0];
rz(0.8008147) q[0];
x q[1];
rz(-1.4919733) q[2];
sx q[2];
rz(-2.6597201) q[2];
sx q[2];
rz(0.80480591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7501523) q[1];
sx q[1];
rz(-1.2160157) q[1];
sx q[1];
rz(-1.9998235) q[1];
x q[2];
rz(-1.1627694) q[3];
sx q[3];
rz(-0.84822908) q[3];
sx q[3];
rz(-0.17150294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1839972) q[2];
sx q[2];
rz(-0.60773578) q[2];
sx q[2];
rz(-0.11873928) q[2];
rz(2.8710098) q[3];
sx q[3];
rz(-2.1515473) q[3];
sx q[3];
rz(-2.9873007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1110693) q[0];
sx q[0];
rz(-1.9707762) q[0];
sx q[0];
rz(2.6861374) q[0];
rz(-2.9083374) q[1];
sx q[1];
rz(-0.74359727) q[1];
sx q[1];
rz(-3.1366248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7457897) q[0];
sx q[0];
rz(-3.0533049) q[0];
sx q[0];
rz(0.39678617) q[0];
rz(1.8894779) q[2];
sx q[2];
rz(-1.2010842) q[2];
sx q[2];
rz(-1.8869417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3525654) q[1];
sx q[1];
rz(-1.066117) q[1];
sx q[1];
rz(2.6093409) q[1];
rz(-1.6065323) q[3];
sx q[3];
rz(-2.7155345) q[3];
sx q[3];
rz(-2.1781875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3576971) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(0.89007393) q[2];
rz(0.87738532) q[3];
sx q[3];
rz(-2.2062517) q[3];
sx q[3];
rz(-2.4963511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18168618) q[0];
sx q[0];
rz(-3.063995) q[0];
sx q[0];
rz(2.087387) q[0];
rz(-0.28610817) q[1];
sx q[1];
rz(-1.0461297) q[1];
sx q[1];
rz(1.2623164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842997) q[0];
sx q[0];
rz(-1.2861006) q[0];
sx q[0];
rz(0.64356743) q[0];
rz(-pi) q[1];
rz(1.8334062) q[2];
sx q[2];
rz(-1.76909) q[2];
sx q[2];
rz(1.8526371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8827829) q[1];
sx q[1];
rz(-1.7534326) q[1];
sx q[1];
rz(2.6983077) q[1];
rz(-0.083857337) q[3];
sx q[3];
rz(-1.0174804) q[3];
sx q[3];
rz(-0.38804873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1288278) q[2];
sx q[2];
rz(-0.78170693) q[2];
sx q[2];
rz(2.8774366) q[2];
rz(2.9014897) q[3];
sx q[3];
rz(-2.0930347) q[3];
sx q[3];
rz(-0.31606328) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964597) q[0];
sx q[0];
rz(-0.75440732) q[0];
sx q[0];
rz(-0.9040133) q[0];
rz(-0.28884197) q[1];
sx q[1];
rz(-1.387351) q[1];
sx q[1];
rz(3.0659061) q[1];
rz(0.31927989) q[2];
sx q[2];
rz(-1.8055331) q[2];
sx q[2];
rz(-1.8100678) q[2];
rz(-0.51863019) q[3];
sx q[3];
rz(-0.082280686) q[3];
sx q[3];
rz(-0.12456457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
