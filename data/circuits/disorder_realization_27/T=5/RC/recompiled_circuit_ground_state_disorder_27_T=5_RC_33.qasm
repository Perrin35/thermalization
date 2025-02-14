OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(-2.5728777) q[0];
rz(1.2110127) q[1];
sx q[1];
rz(-2.14415) q[1];
sx q[1];
rz(-0.2675736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8519605) q[0];
sx q[0];
rz(-1.3697764) q[0];
sx q[0];
rz(-2.5876849) q[0];
x q[1];
rz(1.8873147) q[2];
sx q[2];
rz(-2.7915349) q[2];
sx q[2];
rz(-1.8560662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1700344) q[1];
sx q[1];
rz(-1.5791248) q[1];
sx q[1];
rz(-0.002408601) q[1];
x q[2];
rz(-0.062248793) q[3];
sx q[3];
rz(-1.4302974) q[3];
sx q[3];
rz(1.2195017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25207511) q[2];
sx q[2];
rz(-2.0740261) q[2];
sx q[2];
rz(-1.1313103) q[2];
rz(2.6913397) q[3];
sx q[3];
rz(-0.69142747) q[3];
sx q[3];
rz(-2.1972307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(1.7339647) q[0];
rz(-0.44250008) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(2.5462467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93564088) q[0];
sx q[0];
rz(-1.8186339) q[0];
sx q[0];
rz(-1.3252844) q[0];
rz(-2.6518455) q[2];
sx q[2];
rz(-2.1648266) q[2];
sx q[2];
rz(-0.56383946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41619424) q[1];
sx q[1];
rz(-0.42896118) q[1];
sx q[1];
rz(-3.129175) q[1];
rz(1.1717623) q[3];
sx q[3];
rz(-2.1431987) q[3];
sx q[3];
rz(-2.8295598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2596316) q[2];
sx q[2];
rz(-2.5233848) q[2];
sx q[2];
rz(-2.372443) q[2];
rz(-1.361557) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(-0.59534016) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.896987) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(0.33682522) q[0];
rz(-1.4312076) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(0.097188458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48477706) q[0];
sx q[0];
rz(-2.6069399) q[0];
sx q[0];
rz(-1.1645681) q[0];
x q[1];
rz(1.4732185) q[2];
sx q[2];
rz(-1.9505902) q[2];
sx q[2];
rz(1.9623836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9758476) q[1];
sx q[1];
rz(-0.37384181) q[1];
sx q[1];
rz(-2.7553619) q[1];
rz(-pi) q[2];
rz(-0.52690701) q[3];
sx q[3];
rz(-0.95427536) q[3];
sx q[3];
rz(-1.34672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1691957) q[2];
sx q[2];
rz(-1.7652067) q[2];
sx q[2];
rz(0.81673679) q[2];
rz(-0.9225325) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(1.9492662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794491) q[0];
sx q[0];
rz(-1.8035996) q[0];
sx q[0];
rz(0.86679593) q[0];
rz(1.2358933) q[1];
sx q[1];
rz(-2.0265323) q[1];
sx q[1];
rz(-2.8996276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.73409) q[0];
sx q[0];
rz(-0.76341141) q[0];
sx q[0];
rz(1.9620738) q[0];
rz(-pi) q[1];
rz(1.1461805) q[2];
sx q[2];
rz(-0.17689366) q[2];
sx q[2];
rz(1.572027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82034012) q[1];
sx q[1];
rz(-2.7124321) q[1];
sx q[1];
rz(-0.33351516) q[1];
rz(-0.091598467) q[3];
sx q[3];
rz(-2.6428416) q[3];
sx q[3];
rz(2.2745511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26291004) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(-2.6061457) q[2];
rz(0.49324909) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(-2.6935327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0754452) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(2.2976663) q[0];
rz(1.6663724) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(2.7484238) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1895904) q[0];
sx q[0];
rz(-1.31685) q[0];
sx q[0];
rz(-2.3334731) q[0];
x q[1];
rz(2.8507502) q[2];
sx q[2];
rz(-2.484349) q[2];
sx q[2];
rz(0.90279451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6715309) q[1];
sx q[1];
rz(-1.1137149) q[1];
sx q[1];
rz(1.2600598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0601021) q[3];
sx q[3];
rz(-2.296128) q[3];
sx q[3];
rz(0.368834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4142485) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(-2.7434529) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-2.3530493) q[3];
sx q[3];
rz(1.8383693) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(2.2858802) q[0];
rz(2.4489467) q[1];
sx q[1];
rz(-2.1470862) q[1];
sx q[1];
rz(-2.8725502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97996431) q[0];
sx q[0];
rz(-1.3792896) q[0];
sx q[0];
rz(-3.1167517) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6330209) q[2];
sx q[2];
rz(-2.6341558) q[2];
sx q[2];
rz(-2.1560046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0647886) q[1];
sx q[1];
rz(-0.85569438) q[1];
sx q[1];
rz(0.24311693) q[1];
rz(-pi) q[2];
rz(0.16420047) q[3];
sx q[3];
rz(-1.3815666) q[3];
sx q[3];
rz(-2.7606719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1143703) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(3.031292) q[2];
rz(2.2731764) q[3];
sx q[3];
rz(-1.1766368) q[3];
sx q[3];
rz(1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39349839) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(2.2802343) q[0];
rz(-1.512108) q[1];
sx q[1];
rz(-1.4605099) q[1];
sx q[1];
rz(-2.3177564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0938823) q[0];
sx q[0];
rz(-1.4589696) q[0];
sx q[0];
rz(0.98112962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3917213) q[2];
sx q[2];
rz(-1.2629384) q[2];
sx q[2];
rz(1.4726382) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41125248) q[1];
sx q[1];
rz(-1.3578051) q[1];
sx q[1];
rz(-0.61720444) q[1];
rz(-pi) q[2];
rz(-2.4895913) q[3];
sx q[3];
rz(-1.8240415) q[3];
sx q[3];
rz(0.30919231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5065826) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(3.1363764) q[3];
sx q[3];
rz(-2.3514533) q[3];
sx q[3];
rz(1.691157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.1714627) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(0.18950732) q[0];
rz(-0.7729404) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(0.596284) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8533289) q[0];
sx q[0];
rz(-1.402308) q[0];
sx q[0];
rz(2.4503319) q[0];
rz(-pi) q[1];
rz(-0.79769602) q[2];
sx q[2];
rz(-2.0828649) q[2];
sx q[2];
rz(1.2810117) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.90702) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(1.6638882) q[1];
rz(-pi) q[2];
rz(-0.85875384) q[3];
sx q[3];
rz(-1.158445) q[3];
sx q[3];
rz(1.841973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4160055) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(-1.2131946) q[2];
rz(-2.1554135) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(-1.8825611) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0202494) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(-1.1100618) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(-1.8448255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27889869) q[0];
sx q[0];
rz(-1.4771013) q[0];
sx q[0];
rz(0.91178943) q[0];
rz(1.4234366) q[2];
sx q[2];
rz(-2.7412716) q[2];
sx q[2];
rz(1.6860698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1549266) q[1];
sx q[1];
rz(-2.6441751) q[1];
sx q[1];
rz(-3.0104396) q[1];
rz(-pi) q[2];
rz(-0.66547439) q[3];
sx q[3];
rz(-1.7413119) q[3];
sx q[3];
rz(-2.1297034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0686331) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(-3.0375321) q[2];
rz(-2.0067298) q[3];
sx q[3];
rz(-1.7770504) q[3];
sx q[3];
rz(-2.9141736) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65570152) q[0];
sx q[0];
rz(-2.1567397) q[0];
sx q[0];
rz(0.73053288) q[0];
rz(-2.5841374) q[1];
sx q[1];
rz(-1.971222) q[1];
sx q[1];
rz(2.7117859) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9011079) q[0];
sx q[0];
rz(-2.3798124) q[0];
sx q[0];
rz(-3.1408589) q[0];
x q[1];
rz(-2.5556106) q[2];
sx q[2];
rz(-0.44011099) q[2];
sx q[2];
rz(0.42428478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.046172) q[1];
sx q[1];
rz(-2.7616427) q[1];
sx q[1];
rz(-1.7871961) q[1];
x q[2];
rz(-1.2602706) q[3];
sx q[3];
rz(-0.38703296) q[3];
sx q[3];
rz(-0.68614764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80959117) q[2];
sx q[2];
rz(-2.1198699) q[2];
sx q[2];
rz(0.29279718) q[2];
rz(-3.0012567) q[3];
sx q[3];
rz(-1.0293181) q[3];
sx q[3];
rz(-0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705227) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(0.71939214) q[1];
sx q[1];
rz(-1.2750625) q[1];
sx q[1];
rz(0.19663179) q[1];
rz(1.0631845) q[2];
sx q[2];
rz(-1.4127991) q[2];
sx q[2];
rz(0.71411919) q[2];
rz(3.0075913) q[3];
sx q[3];
rz(-1.8905427) q[3];
sx q[3];
rz(-0.34656634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
