OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(4.1252131) q[1];
sx q[1];
rz(10.317378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8164506) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(1.0613326) q[0];
rz(-2.4721488) q[2];
sx q[2];
rz(-1.1602243) q[2];
sx q[2];
rz(-2.4759811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(2.288504) q[1];
rz(-2.3905972) q[3];
sx q[3];
rz(-1.5392116) q[3];
sx q[3];
rz(-2.3834474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.974568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.5123899) q[0];
sx q[0];
rz(-1.8286684) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8950047) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.8794683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6845219) q[1];
sx q[1];
rz(-1.5091389) q[1];
sx q[1];
rz(2.3350299) q[1];
x q[2];
rz(-2.5490709) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(0.034686397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(-3.068148) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-1.0148369) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944511) q[0];
sx q[0];
rz(-1.6269636) q[0];
sx q[0];
rz(2.4980314) q[0];
x q[1];
rz(-2.0929167) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(0.93333474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2753607) q[1];
sx q[1];
rz(-2.1664157) q[1];
sx q[1];
rz(0.097252107) q[1];
rz(-pi) q[2];
rz(-0.88481836) q[3];
sx q[3];
rz(-1.3664477) q[3];
sx q[3];
rz(-2.200978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(2.5057709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5084002) q[0];
sx q[0];
rz(-2.0661372) q[0];
sx q[0];
rz(0.20423996) q[0];
x q[1];
rz(-1.7556778) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(0.76279574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35866657) q[1];
sx q[1];
rz(-1.7554605) q[1];
sx q[1];
rz(2.3653046) q[1];
rz(-pi) q[2];
rz(-2.8549529) q[3];
sx q[3];
rz(-1.3661824) q[3];
sx q[3];
rz(2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9277966) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(-0.27854663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51019788) q[0];
sx q[0];
rz(-0.83980951) q[0];
sx q[0];
rz(2.1561119) q[0];
rz(2.9859221) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.2200668) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4651471) q[1];
sx q[1];
rz(-0.74063456) q[1];
sx q[1];
rz(-1.0100823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99028011) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(-0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-3.0335398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6022588) q[0];
sx q[0];
rz(-2.5784011) q[0];
sx q[0];
rz(2.5186033) q[0];
rz(-0.33424218) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(-0.7451171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3012645) q[1];
sx q[1];
rz(-0.64126188) q[1];
sx q[1];
rz(-2.458359) q[1];
rz(-1.990854) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(-0.89362842) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-2.4126841) q[0];
sx q[0];
rz(2.6364987) q[0];
rz(1.0412752) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(-1.668001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6210729) q[1];
sx q[1];
rz(-2.2398661) q[1];
sx q[1];
rz(-1.8280162) q[1];
x q[2];
rz(-0.75974792) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(-2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(-2.0513127) q[0];
rz(-3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(0.056218938) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7190785) q[2];
sx q[2];
rz(-2.1170756) q[2];
sx q[2];
rz(0.15020457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1057508) q[1];
sx q[1];
rz(-2.3301947) q[1];
sx q[1];
rz(1.4856542) q[1];
rz(-pi) q[2];
rz(-2.7525547) q[3];
sx q[3];
rz(-0.93718796) q[3];
sx q[3];
rz(0.47215677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(-1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4719452) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(1.9452842) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8413999) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(-1.1635309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3635892) q[1];
sx q[1];
rz(-2.6215141) q[1];
sx q[1];
rz(0.37104721) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2421078) q[3];
sx q[3];
rz(-1.4624819) q[3];
sx q[3];
rz(1.2182106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(2.6616667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23586789) q[0];
sx q[0];
rz(-2.1972482) q[0];
sx q[0];
rz(0.41525526) q[0];
x q[1];
rz(-0.23649044) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(2.1361534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-2.5017782) q[1];
sx q[1];
rz(-2.5261643) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89703538) q[3];
sx q[3];
rz(-1.9978791) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823572) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(0.97933979) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(-1.7469035) q[2];
sx q[2];
rz(-1.0955784) q[2];
sx q[2];
rz(-0.57217862) q[2];
rz(1.3027719) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
