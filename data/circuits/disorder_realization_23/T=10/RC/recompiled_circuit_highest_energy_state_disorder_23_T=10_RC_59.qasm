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
rz(0.5956369) q[0];
sx q[0];
rz(3.5516153) q[0];
sx q[0];
rz(8.6300996) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(5.3310634) q[1];
sx q[1];
rz(7.4330243) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2993752) q[0];
sx q[0];
rz(-1.8646075) q[0];
sx q[0];
rz(3.1042645) q[0];
x q[1];
rz(-0.29157421) q[2];
sx q[2];
rz(-2.2842424) q[2];
sx q[2];
rz(-0.50285646) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1049624) q[1];
sx q[1];
rz(-1.7425907) q[1];
sx q[1];
rz(0.40375945) q[1];
rz(-pi) q[2];
rz(1.8273699) q[3];
sx q[3];
rz(-1.1169516) q[3];
sx q[3];
rz(-0.75405771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6940234) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(2.0726974) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.4890891) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(0.51134837) q[0];
rz(-1.6343575) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(1.7587597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16973142) q[0];
sx q[0];
rz(-1.7771374) q[0];
sx q[0];
rz(1.1792762) q[0];
x q[1];
rz(-2.2425425) q[2];
sx q[2];
rz(-1.8457638) q[2];
sx q[2];
rz(3.0877047) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.4310337) q[1];
sx q[1];
rz(-1.1712499) q[1];
sx q[1];
rz(-2.8053158) q[1];
rz(-pi) q[2];
rz(2.9622127) q[3];
sx q[3];
rz(-0.53032833) q[3];
sx q[3];
rz(1.3190003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4649268) q[2];
sx q[2];
rz(-1.2049462) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(1.0627221) q[3];
sx q[3];
rz(-1.9461742) q[3];
sx q[3];
rz(3.0265871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17387834) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(-0.1097196) q[0];
rz(2.8166855) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(-2.6085764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.479411) q[0];
sx q[0];
rz(-1.4328188) q[0];
sx q[0];
rz(-2.8421287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7092878) q[2];
sx q[2];
rz(-3.084331) q[2];
sx q[2];
rz(-2.6309516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91777201) q[1];
sx q[1];
rz(-0.9448565) q[1];
sx q[1];
rz(-0.84898941) q[1];
rz(-pi) q[2];
rz(2.7044933) q[3];
sx q[3];
rz(-2.282939) q[3];
sx q[3];
rz(2.8959078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.816232) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(1.4556966) q[2];
rz(3.0560737) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(-2.1948309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(-1.5708615) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-3.0991203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149276) q[0];
sx q[0];
rz(-2.6579497) q[0];
sx q[0];
rz(2.7381104) q[0];
rz(-pi) q[1];
rz(-1.4804777) q[2];
sx q[2];
rz(-1.5549763) q[2];
sx q[2];
rz(-2.9243369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3610246) q[1];
sx q[1];
rz(-0.79611049) q[1];
sx q[1];
rz(-1.8382422) q[1];
rz(-1.8117531) q[3];
sx q[3];
rz(-0.78380871) q[3];
sx q[3];
rz(-0.93894053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2461569) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(-0.79143381) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(-2.4401276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0303845) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(2.9682888) q[0];
rz(-0.93470848) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(2.8876143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58502176) q[0];
sx q[0];
rz(-2.3818172) q[0];
sx q[0];
rz(1.7346481) q[0];
rz(-2.6298275) q[2];
sx q[2];
rz(-2.4825271) q[2];
sx q[2];
rz(2.2816254) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3601927) q[1];
sx q[1];
rz(-1.5500229) q[1];
sx q[1];
rz(-2.4482577) q[1];
rz(-1.41296) q[3];
sx q[3];
rz(-0.39208347) q[3];
sx q[3];
rz(-0.81721701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0678593) q[2];
sx q[2];
rz(-0.24418712) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(-0.67617792) q[3];
sx q[3];
rz(-1.9192326) q[3];
sx q[3];
rz(0.079782709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22236958) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(-2.4134912) q[0];
rz(1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(2.5808835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8122305) q[0];
sx q[0];
rz(-1.6971611) q[0];
sx q[0];
rz(1.6362599) q[0];
x q[1];
rz(2.9496983) q[2];
sx q[2];
rz(-2.5835496) q[2];
sx q[2];
rz(3.0241338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.079472629) q[1];
sx q[1];
rz(-1.9626856) q[1];
sx q[1];
rz(-2.9074039) q[1];
rz(-pi) q[2];
rz(0.57413103) q[3];
sx q[3];
rz(-3.0256292) q[3];
sx q[3];
rz(2.9812893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59549436) q[2];
sx q[2];
rz(-1.6807154) q[2];
sx q[2];
rz(-2.0687436) q[2];
rz(-2.8192375) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(0.391092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748465) q[0];
sx q[0];
rz(-1.5006737) q[0];
sx q[0];
rz(2.7333976) q[0];
rz(0.32632581) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(1.6654642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832907) q[0];
sx q[0];
rz(-1.5314449) q[0];
sx q[0];
rz(1.395306) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23952837) q[2];
sx q[2];
rz(-2.1210741) q[2];
sx q[2];
rz(0.36877647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5567814) q[1];
sx q[1];
rz(-1.4540392) q[1];
sx q[1];
rz(-0.054971855) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9093698) q[3];
sx q[3];
rz(-1.4055863) q[3];
sx q[3];
rz(-2.8286407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8866715) q[2];
sx q[2];
rz(-1.485606) q[2];
sx q[2];
rz(-2.2754748) q[2];
rz(-2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(2.2542663) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0199468) q[0];
sx q[0];
rz(-0.045504657) q[0];
sx q[0];
rz(1.2787) q[0];
rz(-1.0026898) q[1];
sx q[1];
rz(-1.2607144) q[1];
sx q[1];
rz(0.444828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0507752) q[0];
sx q[0];
rz(-1.7184166) q[0];
sx q[0];
rz(2.0992797) q[0];
x q[1];
rz(2.41909) q[2];
sx q[2];
rz(-1.2151657) q[2];
sx q[2];
rz(0.96273811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.63193601) q[1];
sx q[1];
rz(-0.77101123) q[1];
sx q[1];
rz(-1.7421175) q[1];
rz(-1.3108389) q[3];
sx q[3];
rz(-1.9366067) q[3];
sx q[3];
rz(-1.8577675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7030846) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-2.6355696) q[2];
rz(2.1972726) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.849702) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(-3.131409) q[0];
rz(-0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-0.082286509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1707538) q[0];
sx q[0];
rz(-1.0910831) q[0];
sx q[0];
rz(-1.1520349) q[0];
rz(-pi) q[1];
rz(-0.89110969) q[2];
sx q[2];
rz(-1.3902877) q[2];
sx q[2];
rz(-1.0918822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.391606) q[1];
sx q[1];
rz(-0.22400236) q[1];
sx q[1];
rz(0.27952607) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27307053) q[3];
sx q[3];
rz(-2.4585173) q[3];
sx q[3];
rz(-3.0103252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8102707) q[2];
sx q[2];
rz(-1.6931345) q[2];
sx q[2];
rz(-2.4324379) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(2.6813996) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3191147) q[0];
sx q[0];
rz(-0.8304441) q[0];
sx q[0];
rz(0.6293695) q[0];
rz(1.7259701) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.9633044) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30692682) q[0];
sx q[0];
rz(-0.7531868) q[0];
sx q[0];
rz(-0.96145328) q[0];
rz(-pi) q[1];
rz(-1.0683996) q[2];
sx q[2];
rz(-1.2594688) q[2];
sx q[2];
rz(-2.5466998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3523184) q[1];
sx q[1];
rz(-1.6802009) q[1];
sx q[1];
rz(2.3477702) q[1];
rz(-pi) q[2];
rz(-0.72145907) q[3];
sx q[3];
rz(-1.1260403) q[3];
sx q[3];
rz(-0.12192425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(2.3229522) q[2];
rz(1.5107752) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(-2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4928987) q[0];
sx q[0];
rz(-1.10981) q[0];
sx q[0];
rz(1.4640402) q[0];
rz(1.875444) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(-1.2513592) q[2];
sx q[2];
rz(-0.70079679) q[2];
sx q[2];
rz(0.23637017) q[2];
rz(2.6880445) q[3];
sx q[3];
rz(-0.20487479) q[3];
sx q[3];
rz(0.80036001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
