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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(-2.3469143) q[0];
rz(2.088264) q[1];
sx q[1];
rz(-2.1894708) q[1];
sx q[1];
rz(1.9917537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7139706) q[0];
sx q[0];
rz(-2.8454878) q[0];
sx q[0];
rz(1.6935191) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29157421) q[2];
sx q[2];
rz(-0.85735029) q[2];
sx q[2];
rz(0.50285646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6803327) q[1];
sx q[1];
rz(-1.1733175) q[1];
sx q[1];
rz(-1.3843126) q[1];
rz(2.6618979) q[3];
sx q[3];
rz(-0.51691662) q[3];
sx q[3];
rz(1.848288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4475693) q[2];
sx q[2];
rz(-2.1489216) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(2.0726974) q[3];
sx q[3];
rz(-1.5690208) q[3];
sx q[3];
rz(1.6525035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25817961) q[0];
sx q[0];
rz(-3.1386107) q[0];
sx q[0];
rz(2.6302443) q[0];
rz(-1.5072352) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(-1.3828329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718612) q[0];
sx q[0];
rz(-1.3644553) q[0];
sx q[0];
rz(-1.9623164) q[0];
rz(2.2425425) q[2];
sx q[2];
rz(-1.2958289) q[2];
sx q[2];
rz(-0.053887966) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0046151) q[1];
sx q[1];
rz(-1.8796693) q[1];
sx q[1];
rz(-1.9914133) q[1];
rz(-pi) q[2];
rz(-1.6750349) q[3];
sx q[3];
rz(-1.0498507) q[3];
sx q[3];
rz(-1.5261861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4649268) q[2];
sx q[2];
rz(-1.2049462) q[2];
sx q[2];
rz(-3.0386772) q[2];
rz(-1.0627221) q[3];
sx q[3];
rz(-1.9461742) q[3];
sx q[3];
rz(-3.0265871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17387834) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(-3.0318731) q[0];
rz(-2.8166855) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(-0.5330162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1905381) q[0];
sx q[0];
rz(-1.2742654) q[0];
sx q[0];
rz(1.426479) q[0];
rz(-1.5948086) q[2];
sx q[2];
rz(-1.5188076) q[2];
sx q[2];
rz(2.1980224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2238206) q[1];
sx q[1];
rz(-2.1967362) q[1];
sx q[1];
rz(2.2926032) q[1];
rz(2.332052) q[3];
sx q[3];
rz(-1.8969826) q[3];
sx q[3];
rz(1.0287788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32536062) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(1.4556966) q[2];
rz(-0.085518941) q[3];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(1.5708615) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-0.042472366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6149276) q[0];
sx q[0];
rz(-0.4836429) q[0];
sx q[0];
rz(-0.40348224) q[0];
x q[1];
rz(-1.3971523) q[2];
sx q[2];
rz(-0.091689907) q[2];
sx q[2];
rz(1.9609811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5407357) q[1];
sx q[1];
rz(-1.3807978) q[1];
sx q[1];
rz(-2.3488087) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23353429) q[3];
sx q[3];
rz(-2.326205) q[3];
sx q[3];
rz(2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.3763206) q[2];
sx q[2];
rz(2.8975272) q[2];
rz(-0.79143381) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1112082) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(-2.9682888) q[0];
rz(-0.93470848) q[1];
sx q[1];
rz(-2.2917031) q[1];
sx q[1];
rz(-0.25397837) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365216) q[0];
sx q[0];
rz(-1.6833841) q[0];
sx q[0];
rz(-2.3238411) q[0];
rz(-0.59400174) q[2];
sx q[2];
rz(-1.8753759) q[2];
sx q[2];
rz(2.0128743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22786247) q[1];
sx q[1];
rz(-0.87764064) q[1];
sx q[1];
rz(-1.5437897) q[1];
rz(1.41296) q[3];
sx q[3];
rz(-2.7495092) q[3];
sx q[3];
rz(-0.81721701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0678593) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(1.4402639) q[2];
rz(-0.67617792) q[3];
sx q[3];
rz(-1.22236) q[3];
sx q[3];
rz(3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9192231) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(0.72810143) q[0];
rz(-1.1657731) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(-2.5808835) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3293622) q[0];
sx q[0];
rz(-1.4444315) q[0];
sx q[0];
rz(1.6362599) q[0];
rz(-pi) q[1];
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
rz(0.079472629) q[1];
sx q[1];
rz(-1.9626856) q[1];
sx q[1];
rz(2.9074039) q[1];
rz(-3.0440935) q[3];
sx q[3];
rz(-1.6336771) q[3];
sx q[3];
rz(1.9815552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5460983) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(-1.072849) q[2];
rz(0.32235518) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(-2.7505007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.40819502) q[0];
rz(-0.32632581) q[1];
sx q[1];
rz(-0.34714454) q[1];
sx q[1];
rz(1.6654642) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30551874) q[0];
sx q[0];
rz(-1.7461494) q[0];
sx q[0];
rz(-0.039964635) q[0];
rz(-pi) q[1];
rz(-1.0075241) q[2];
sx q[2];
rz(-1.7744229) q[2];
sx q[2];
rz(1.8125535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.020425107) q[1];
sx q[1];
rz(-1.6253935) q[1];
sx q[1];
rz(1.6877285) q[1];
rz(-pi) q[2];
rz(2.9093698) q[3];
sx q[3];
rz(-1.4055863) q[3];
sx q[3];
rz(-2.8286407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8866715) q[2];
sx q[2];
rz(-1.485606) q[2];
sx q[2];
rz(2.2754748) q[2];
rz(0.73650375) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(2.2542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0199468) q[0];
sx q[0];
rz(-0.045504657) q[0];
sx q[0];
rz(1.2787) q[0];
rz(1.0026898) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-2.6967646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0908175) q[0];
sx q[0];
rz(-1.4231761) q[0];
sx q[0];
rz(-2.0992797) q[0];
x q[1];
rz(1.1110531) q[2];
sx q[2];
rz(-2.2393949) q[2];
sx q[2];
rz(2.8313314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0792744) q[1];
sx q[1];
rz(-1.6898815) q[1];
sx q[1];
rz(2.3344385) q[1];
rz(-pi) q[2];
rz(0.37737198) q[3];
sx q[3];
rz(-1.8131953) q[3];
sx q[3];
rz(-2.9494772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-0.50602305) q[2];
rz(0.94432008) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(2.4054312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29189062) q[0];
sx q[0];
rz(-2.1493752) q[0];
sx q[0];
rz(0.010183656) q[0];
rz(0.61095515) q[1];
sx q[1];
rz(-0.98038951) q[1];
sx q[1];
rz(-0.082286509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9708389) q[0];
sx q[0];
rz(-2.0505095) q[0];
sx q[0];
rz(1.1520349) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23046979) q[2];
sx q[2];
rz(-2.2374212) q[2];
sx q[2];
rz(2.8067592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4520784) q[1];
sx q[1];
rz(-1.6321215) q[1];
sx q[1];
rz(2.9260127) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8685221) q[3];
sx q[3];
rz(-0.6830754) q[3];
sx q[3];
rz(0.13126743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8102707) q[2];
sx q[2];
rz(-1.6931345) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-0.63153657) q[3];
sx q[3];
rz(-2.6813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.822478) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(-1.7259701) q[1];
sx q[1];
rz(-1.4868088) q[1];
sx q[1];
rz(1.1782882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45635763) q[0];
sx q[0];
rz(-2.166232) q[0];
sx q[0];
rz(-2.6490982) q[0];
rz(-pi) q[1];
rz(-2.1599102) q[2];
sx q[2];
rz(-0.58393634) q[2];
sx q[2];
rz(-2.6743025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32543031) q[1];
sx q[1];
rz(-0.79968444) q[1];
sx q[1];
rz(2.9887448) q[1];
x q[2];
rz(0.72145907) q[3];
sx q[3];
rz(-1.1260403) q[3];
sx q[3];
rz(0.12192425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87469953) q[2];
sx q[2];
rz(-1.5692254) q[2];
sx q[2];
rz(0.81864041) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4928987) q[0];
sx q[0];
rz(-2.0317827) q[0];
sx q[0];
rz(-1.6775525) q[0];
rz(1.2661487) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(-1.8902334) q[2];
sx q[2];
rz(-2.4407959) q[2];
sx q[2];
rz(-2.9052225) q[2];
rz(-0.18465445) q[3];
sx q[3];
rz(-1.481537) q[3];
sx q[3];
rz(2.8164345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
