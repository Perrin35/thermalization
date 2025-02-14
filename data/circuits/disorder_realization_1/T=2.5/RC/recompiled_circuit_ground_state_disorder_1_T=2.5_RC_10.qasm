OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(-2.1362526) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(-2.7856196) q[1];
sx q[1];
rz(2.2495143) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1482502) q[0];
sx q[0];
rz(-0.82776946) q[0];
sx q[0];
rz(-1.7048852) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96718915) q[2];
sx q[2];
rz(-1.4481067) q[2];
sx q[2];
rz(-1.762378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6309384) q[1];
sx q[1];
rz(-1.0882241) q[1];
sx q[1];
rz(1.6810745) q[1];
rz(-pi) q[2];
rz(-2.0964811) q[3];
sx q[3];
rz(-1.0307923) q[3];
sx q[3];
rz(-1.857855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4165108) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(1.4408646) q[2];
rz(0.95669389) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956534) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(2.0126427) q[0];
rz(-2.8961862) q[1];
sx q[1];
rz(-1.8281728) q[1];
sx q[1];
rz(-2.479877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8186057) q[0];
sx q[0];
rz(-1.855369) q[0];
sx q[0];
rz(-0.63708441) q[0];
rz(-pi) q[1];
rz(-2.6241669) q[2];
sx q[2];
rz(-2.708507) q[2];
sx q[2];
rz(-0.171207) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0894818) q[1];
sx q[1];
rz(-2.8683337) q[1];
sx q[1];
rz(0.11788003) q[1];
rz(-pi) q[2];
rz(-1.8850967) q[3];
sx q[3];
rz(-2.8366025) q[3];
sx q[3];
rz(1.5105607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5976065) q[2];
sx q[2];
rz(-2.643955) q[2];
sx q[2];
rz(-2.0987161) q[2];
rz(1.6889702) q[3];
sx q[3];
rz(-2.2904604) q[3];
sx q[3];
rz(1.1037306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76688981) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(1.8324628) q[0];
rz(2.6922928) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(-1.4264533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9096722) q[0];
sx q[0];
rz(-1.8350775) q[0];
sx q[0];
rz(2.8876165) q[0];
rz(1.6108247) q[2];
sx q[2];
rz(-0.42280254) q[2];
sx q[2];
rz(-1.1943447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92628463) q[1];
sx q[1];
rz(-2.4325772) q[1];
sx q[1];
rz(3.0224817) q[1];
x q[2];
rz(1.8934694) q[3];
sx q[3];
rz(-2.8827658) q[3];
sx q[3];
rz(-0.62773529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36491498) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(0.066369973) q[2];
rz(-2.5868609) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(-1.6248645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9471112) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(-2.3696005) q[0];
rz(-0.66328612) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(3.0139121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6948815) q[0];
sx q[0];
rz(-2.237473) q[0];
sx q[0];
rz(-1.5645909) q[0];
rz(-2.5721385) q[2];
sx q[2];
rz(-2.5094911) q[2];
sx q[2];
rz(0.60213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3290594) q[1];
sx q[1];
rz(-2.529728) q[1];
sx q[1];
rz(-0.81176035) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7112682) q[3];
sx q[3];
rz(-1.7146304) q[3];
sx q[3];
rz(2.6324684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53106236) q[2];
sx q[2];
rz(-0.9610815) q[2];
sx q[2];
rz(0.5985716) q[2];
rz(-0.5851723) q[3];
sx q[3];
rz(-1.3734615) q[3];
sx q[3];
rz(-0.055805512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35036206) q[0];
sx q[0];
rz(-1.5086011) q[0];
sx q[0];
rz(-2.1685261) q[0];
rz(-2.9452501) q[1];
sx q[1];
rz(-1.1390431) q[1];
sx q[1];
rz(-1.6729209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356883) q[0];
sx q[0];
rz(-2.2668512) q[0];
sx q[0];
rz(2.619834) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0650915) q[2];
sx q[2];
rz(-1.7943766) q[2];
sx q[2];
rz(2.2044942) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1717559) q[1];
sx q[1];
rz(-1.7431269) q[1];
sx q[1];
rz(1.8755765) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0319388) q[3];
sx q[3];
rz(-1.6183637) q[3];
sx q[3];
rz(-2.4470234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4520182) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(0.16239521) q[2];
rz(2.0116122) q[3];
sx q[3];
rz(-2.2584848) q[3];
sx q[3];
rz(-1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3910386) q[0];
sx q[0];
rz(-0.51893187) q[0];
sx q[0];
rz(-0.050405141) q[0];
rz(0.16818908) q[1];
sx q[1];
rz(-1.5610118) q[1];
sx q[1];
rz(-2.2680297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03913232) q[0];
sx q[0];
rz(-0.96193571) q[0];
sx q[0];
rz(1.2723421) q[0];
x q[1];
rz(-1.8040405) q[2];
sx q[2];
rz(-2.9961259) q[2];
sx q[2];
rz(-1.5158397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7616523) q[1];
sx q[1];
rz(-0.72387513) q[1];
sx q[1];
rz(-0.70751247) q[1];
rz(-0.29804067) q[3];
sx q[3];
rz(-1.3750769) q[3];
sx q[3];
rz(0.51575553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1217338) q[2];
sx q[2];
rz(-2.9254318) q[2];
sx q[2];
rz(0.38062322) q[2];
rz(0.56626433) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.5009163) q[0];
sx q[0];
rz(-2.930142) q[0];
sx q[0];
rz(0.40263116) q[0];
rz(-1.7598049) q[1];
sx q[1];
rz(-2.333162) q[1];
sx q[1];
rz(3.0852539) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1884679) q[0];
sx q[0];
rz(-1.2446277) q[0];
sx q[0];
rz(-1.1475546) q[0];
x q[1];
rz(1.4903421) q[2];
sx q[2];
rz(-0.47199303) q[2];
sx q[2];
rz(2.9126963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2353387) q[1];
sx q[1];
rz(-1.5829931) q[1];
sx q[1];
rz(-1.0066628) q[1];
x q[2];
rz(-0.67752892) q[3];
sx q[3];
rz(-1.1544246) q[3];
sx q[3];
rz(-2.8512521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23400433) q[2];
sx q[2];
rz(-1.7981217) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(3.0808926) q[3];
sx q[3];
rz(-1.3541636) q[3];
sx q[3];
rz(-0.48262706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.56977) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(1.7806336) q[0];
rz(0.743615) q[1];
sx q[1];
rz(-1.4185602) q[1];
sx q[1];
rz(0.13882151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1847938) q[0];
sx q[0];
rz(-2.0763872) q[0];
sx q[0];
rz(-2.5661945) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97644225) q[2];
sx q[2];
rz(-1.9658372) q[2];
sx q[2];
rz(-2.3740685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3070345) q[1];
sx q[1];
rz(-1.2456248) q[1];
sx q[1];
rz(2.9284655) q[1];
rz(-pi) q[2];
rz(1.9622063) q[3];
sx q[3];
rz(-1.7633121) q[3];
sx q[3];
rz(0.26993902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(0.71411258) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(-0.97904557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7239083) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-3.0133001) q[0];
rz(2.7543606) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(-1.8181575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0751064) q[0];
sx q[0];
rz(-1.3818701) q[0];
sx q[0];
rz(0.099204258) q[0];
rz(-pi) q[1];
rz(-2.5495731) q[2];
sx q[2];
rz(-0.97416234) q[2];
sx q[2];
rz(-0.48812619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.379017) q[1];
sx q[1];
rz(-1.5847995) q[1];
sx q[1];
rz(-1.714731) q[1];
rz(0.72800379) q[3];
sx q[3];
rz(-1.7204855) q[3];
sx q[3];
rz(2.585175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6748176) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(-2.9368994) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650763) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(-0.21328558) q[0];
rz(-0.22008303) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(1.782104) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3228139) q[0];
sx q[0];
rz(-1.5548163) q[0];
sx q[0];
rz(-0.009601618) q[0];
x q[1];
rz(-0.47021265) q[2];
sx q[2];
rz(-2.6192952) q[2];
sx q[2];
rz(3.0015415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6776646) q[1];
sx q[1];
rz(-1.0654133) q[1];
sx q[1];
rz(-2.7457854) q[1];
x q[2];
rz(2.241019) q[3];
sx q[3];
rz(-1.9346721) q[3];
sx q[3];
rz(-2.8577141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.517211) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(-0.33779302) q[2];
rz(-1.7289915) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74078858) q[0];
sx q[0];
rz(-2.3176226) q[0];
sx q[0];
rz(1.2299706) q[0];
rz(0.38372718) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(-2.7376851) q[2];
sx q[2];
rz(-1.2474893) q[2];
sx q[2];
rz(2.9328109) q[2];
rz(2.0819584) q[3];
sx q[3];
rz(-1.3564902) q[3];
sx q[3];
rz(-1.0839957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
