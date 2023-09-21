OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6807251) q[0];
sx q[0];
rz(-1.4493677) q[0];
sx q[0];
rz(-0.073343883) q[0];
x q[1];
rz(3.112552) q[2];
sx q[2];
rz(-1.8543058) q[2];
sx q[2];
rz(-3.0344506) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.764896) q[1];
sx q[1];
rz(-2.2590859) q[1];
sx q[1];
rz(-2.1535758) q[1];
rz(-pi) q[2];
x q[2];
rz(1.430106) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(-1.5593485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.010014023) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-0.28796089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44293091) q[0];
sx q[0];
rz(-1.3804111) q[0];
sx q[0];
rz(-0.18106826) q[0];
rz(2.1638853) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(-1.5175982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49111734) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(2.526545) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35034758) q[3];
sx q[3];
rz(-1.6068034) q[3];
sx q[3];
rz(-2.2909174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59445375) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(0.95169383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0591473) q[0];
sx q[0];
rz(-1.379231) q[0];
sx q[0];
rz(-0.092415718) q[0];
rz(-pi) q[1];
rz(-1.9556324) q[2];
sx q[2];
rz(-0.54883146) q[2];
sx q[2];
rz(1.8897111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23783319) q[1];
sx q[1];
rz(-2.2888498) q[1];
sx q[1];
rz(-0.41463931) q[1];
x q[2];
rz(-2.6056616) q[3];
sx q[3];
rz(-1.0169573) q[3];
sx q[3];
rz(-1.904408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13016985) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(2.7974131) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(0.56604958) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.4973649) q[0];
sx q[0];
rz(1.244506) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(2.7688162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5107721) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(-2.3833016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1431085) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(2.1172303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43416926) q[1];
sx q[1];
rz(-1.8137099) q[1];
sx q[1];
rz(2.0704107) q[1];
rz(-pi) q[2];
rz(2.8126206) q[3];
sx q[3];
rz(-1.1165459) q[3];
sx q[3];
rz(0.3262375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(2.2367031) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.6624517) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(-1.6632535) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6144845) q[0];
sx q[0];
rz(-1.2687578) q[0];
sx q[0];
rz(-1.8943647) q[0];
x q[1];
rz(0.34110951) q[2];
sx q[2];
rz(-2.038539) q[2];
sx q[2];
rz(1.9175921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6492669) q[1];
sx q[1];
rz(-1.6345134) q[1];
sx q[1];
rz(-1.4953514) q[1];
rz(-1.4388496) q[3];
sx q[3];
rz(-2.2852516) q[3];
sx q[3];
rz(0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(0.21970704) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.7664849) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.188376) q[0];
sx q[0];
rz(-1.2326476) q[0];
sx q[0];
rz(-1.2921278) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0639406) q[2];
sx q[2];
rz(-0.75040557) q[2];
sx q[2];
rz(-0.99496182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82735862) q[1];
sx q[1];
rz(-1.595256) q[1];
sx q[1];
rz(2.7315061) q[1];
rz(-pi) q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-2.0241963) q[3];
sx q[3];
rz(0.51844937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(0.19255157) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(0.87316978) q[0];
rz(2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(3.057664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5556363) q[0];
sx q[0];
rz(-1.7935828) q[0];
sx q[0];
rz(1.0046093) q[0];
rz(-pi) q[1];
x q[1];
rz(2.601868) q[2];
sx q[2];
rz(-2.735609) q[2];
sx q[2];
rz(-0.27137953) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9675688) q[1];
sx q[1];
rz(-0.98185437) q[1];
sx q[1];
rz(-2.3322361) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84079068) q[3];
sx q[3];
rz(-1.519324) q[3];
sx q[3];
rz(-1.3703025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(2.7111354) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51171821) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5686533) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(0.22139876) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68265712) q[2];
sx q[2];
rz(-0.97634146) q[2];
sx q[2];
rz(-0.97460954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9981873) q[1];
sx q[1];
rz(-2.7618976) q[1];
sx q[1];
rz(-3.0642964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66611992) q[3];
sx q[3];
rz(-1.564581) q[3];
sx q[3];
rz(-2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-0.60493207) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9873001) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(-2.212021) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(0.5756793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8262779) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(-0.6420352) q[0];
rz(0.11328463) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(1.726113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8262144) q[1];
sx q[1];
rz(-1.4634906) q[1];
sx q[1];
rz(0.025749287) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9541295) q[3];
sx q[3];
rz(-2.1037357) q[3];
sx q[3];
rz(2.370358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(-0.17045505) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-2.5073994) q[0];
rz(0.028907396) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432188) q[0];
sx q[0];
rz(-1.507483) q[0];
sx q[0];
rz(-1.4988585) q[0];
x q[1];
rz(-2.3949941) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.6105086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58418729) q[1];
sx q[1];
rz(-0.69328847) q[1];
sx q[1];
rz(2.8597945) q[1];
x q[2];
rz(1.2877527) q[3];
sx q[3];
rz(-1.9962365) q[3];
sx q[3];
rz(2.2417559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9344899) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(2.2864441) q[2];
sx q[2];
rz(-0.48968857) q[2];
sx q[2];
rz(-0.80077632) q[2];
rz(0.58892693) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];