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
rz(-1.7419185) q[0];
sx q[0];
rz(-1.4631441) q[0];
sx q[0];
rz(1.4611257) q[0];
rz(-0.48179102) q[1];
sx q[1];
rz(2.3269589) q[1];
sx q[1];
rz(8.485312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2964824) q[0];
sx q[0];
rz(-0.044707693) q[0];
sx q[0];
rz(1.787767) q[0];
x q[1];
rz(1.9389802) q[2];
sx q[2];
rz(-1.069931) q[2];
sx q[2];
rz(-0.11239582) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71684696) q[1];
sx q[1];
rz(-1.5813649) q[1];
sx q[1];
rz(-1.6109857) q[1];
rz(-1.788673) q[3];
sx q[3];
rz(-1.9697726) q[3];
sx q[3];
rz(0.50673317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3664832) q[2];
sx q[2];
rz(-1.4668239) q[2];
sx q[2];
rz(0.62466204) q[2];
rz(2.1810253) q[3];
sx q[3];
rz(-1.6865691) q[3];
sx q[3];
rz(-0.87489405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7333154) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(-0.1917924) q[0];
rz(-2.5674112) q[1];
sx q[1];
rz(-1.6185113) q[1];
sx q[1];
rz(-0.40723732) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005368) q[0];
sx q[0];
rz(-1.558702) q[0];
sx q[0];
rz(-1.7465792) q[0];
rz(-pi) q[1];
rz(-1.919433) q[2];
sx q[2];
rz(-2.9538909) q[2];
sx q[2];
rz(-2.3003464) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6169238) q[1];
sx q[1];
rz(-1.7546931) q[1];
sx q[1];
rz(2.6016629) q[1];
x q[2];
rz(-1.7329747) q[3];
sx q[3];
rz(-0.50807488) q[3];
sx q[3];
rz(-1.3401507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.189397) q[2];
sx q[2];
rz(-1.962683) q[2];
sx q[2];
rz(2.4883545) q[2];
rz(0.87853986) q[3];
sx q[3];
rz(-2.0432751) q[3];
sx q[3];
rz(3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.99881309) q[0];
sx q[0];
rz(-1.4343364) q[0];
sx q[0];
rz(2.0913731) q[0];
rz(2.5255919) q[1];
sx q[1];
rz(-1.417825) q[1];
sx q[1];
rz(2.2284257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8693059) q[0];
sx q[0];
rz(-1.4014295) q[0];
sx q[0];
rz(-0.10793751) q[0];
rz(-3.0229125) q[2];
sx q[2];
rz(-1.8556229) q[2];
sx q[2];
rz(0.26155805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29047295) q[1];
sx q[1];
rz(-0.5284068) q[1];
sx q[1];
rz(1.0747133) q[1];
rz(-pi) q[2];
rz(2.2596879) q[3];
sx q[3];
rz(-1.8671452) q[3];
sx q[3];
rz(2.0810082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22127613) q[2];
sx q[2];
rz(-0.90711275) q[2];
sx q[2];
rz(-3.082412) q[2];
rz(-2.3024998) q[3];
sx q[3];
rz(-1.2181166) q[3];
sx q[3];
rz(0.86735094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686104) q[0];
sx q[0];
rz(-2.1324069) q[0];
sx q[0];
rz(2.7918145) q[0];
rz(-1.7971669) q[1];
sx q[1];
rz(-0.58660048) q[1];
sx q[1];
rz(3.0288568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7737425) q[0];
sx q[0];
rz(-2.5597456) q[0];
sx q[0];
rz(0.655158) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2069682) q[2];
sx q[2];
rz(-0.69279659) q[2];
sx q[2];
rz(-0.39510228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4055914) q[1];
sx q[1];
rz(-1.2697769) q[1];
sx q[1];
rz(-2.244414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3808718) q[3];
sx q[3];
rz(-0.84555999) q[3];
sx q[3];
rz(1.4694422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.747644) q[2];
sx q[2];
rz(-1.5851721) q[2];
sx q[2];
rz(-0.057223884) q[2];
rz(-1.8852437) q[3];
sx q[3];
rz(-0.57259721) q[3];
sx q[3];
rz(-0.73067874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2644065) q[0];
sx q[0];
rz(-1.1671966) q[0];
sx q[0];
rz(2.4315244) q[0];
rz(0.26142985) q[1];
sx q[1];
rz(-1.556004) q[1];
sx q[1];
rz(2.1905621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9316021) q[0];
sx q[0];
rz(-1.8977802) q[0];
sx q[0];
rz(-2.6698584) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5102875) q[2];
sx q[2];
rz(-1.1677051) q[2];
sx q[2];
rz(-0.20353488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1360395) q[1];
sx q[1];
rz(-0.73302089) q[1];
sx q[1];
rz(-2.8366778) q[1];
x q[2];
rz(-1.705549) q[3];
sx q[3];
rz(-2.5042078) q[3];
sx q[3];
rz(-2.4911043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49932536) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(-0.30936852) q[2];
rz(-1.0950836) q[3];
sx q[3];
rz(-1.1284004) q[3];
sx q[3];
rz(-2.7373718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9139022) q[0];
sx q[0];
rz(-2.9588283) q[0];
sx q[0];
rz(2.2392654) q[0];
rz(3.0790216) q[1];
sx q[1];
rz(-1.8726655) q[1];
sx q[1];
rz(1.2455469) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093687261) q[0];
sx q[0];
rz(-1.416668) q[0];
sx q[0];
rz(2.4612294) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1104711) q[2];
sx q[2];
rz(-1.6189685) q[2];
sx q[2];
rz(-2.0210514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77038232) q[1];
sx q[1];
rz(-2.4744445) q[1];
sx q[1];
rz(-0.70828153) q[1];
rz(-pi) q[2];
rz(2.4687103) q[3];
sx q[3];
rz(-1.9599293) q[3];
sx q[3];
rz(-1.0997888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41783276) q[2];
sx q[2];
rz(-1.2699026) q[2];
sx q[2];
rz(-1.343824) q[2];
rz(-1.8580681) q[3];
sx q[3];
rz(-2.696974) q[3];
sx q[3];
rz(-0.72492391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.035324) q[0];
sx q[0];
rz(-0.037287354) q[0];
sx q[0];
rz(2.7653747) q[0];
rz(-2.6694934) q[1];
sx q[1];
rz(-0.68196982) q[1];
sx q[1];
rz(0.45904407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58595198) q[0];
sx q[0];
rz(-1.1496854) q[0];
sx q[0];
rz(2.3512406) q[0];
rz(1.3960346) q[2];
sx q[2];
rz(-1.1272001) q[2];
sx q[2];
rz(-0.52973739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7524779) q[1];
sx q[1];
rz(-1.6387617) q[1];
sx q[1];
rz(-0.1631921) q[1];
rz(2.543942) q[3];
sx q[3];
rz(-1.581218) q[3];
sx q[3];
rz(-0.29779336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0692733) q[2];
sx q[2];
rz(-2.4299712) q[2];
sx q[2];
rz(-0.35815987) q[2];
rz(0.88851309) q[3];
sx q[3];
rz(-1.0212967) q[3];
sx q[3];
rz(2.2108938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0446123) q[0];
sx q[0];
rz(-0.29878578) q[0];
sx q[0];
rz(-0.97691798) q[0];
rz(-0.75374976) q[1];
sx q[1];
rz(-1.6981373) q[1];
sx q[1];
rz(1.6389729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328505) q[0];
sx q[0];
rz(-0.70510222) q[0];
sx q[0];
rz(0.1211959) q[0];
rz(-pi) q[1];
rz(2.4719878) q[2];
sx q[2];
rz(-1.1455702) q[2];
sx q[2];
rz(-2.2086672) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9257207) q[1];
sx q[1];
rz(-1.873107) q[1];
sx q[1];
rz(-2.4420101) q[1];
x q[2];
rz(1.8863932) q[3];
sx q[3];
rz(-0.64794174) q[3];
sx q[3];
rz(1.2256988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94540191) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(-0.84575829) q[2];
rz(2.3614007) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(-2.5030524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3188401) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(-2.7986797) q[0];
rz(-2.137939) q[1];
sx q[1];
rz(-1.2316848) q[1];
sx q[1];
rz(2.4249605) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9166853) q[0];
sx q[0];
rz(-1.2222079) q[0];
sx q[0];
rz(2.3681568) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49375294) q[2];
sx q[2];
rz(-1.8519562) q[2];
sx q[2];
rz(-1.2501095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9327859) q[1];
sx q[1];
rz(-2.1541671) q[1];
sx q[1];
rz(1.2482743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9699205) q[3];
sx q[3];
rz(-0.93825996) q[3];
sx q[3];
rz(0.61033953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0922682) q[2];
sx q[2];
rz(-1.5169787) q[2];
sx q[2];
rz(-1.0005023) q[2];
rz(-0.81691027) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(0.36417714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889121) q[0];
sx q[0];
rz(-0.93343821) q[0];
sx q[0];
rz(-0.77891427) q[0];
rz(2.3498416) q[1];
sx q[1];
rz(-1.9154895) q[1];
sx q[1];
rz(-2.7955999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2591922) q[0];
sx q[0];
rz(-2.0171875) q[0];
sx q[0];
rz(-1.6631699) q[0];
rz(-1.6093639) q[2];
sx q[2];
rz(-2.8380605) q[2];
sx q[2];
rz(-0.26798778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3521063) q[1];
sx q[1];
rz(-0.54692422) q[1];
sx q[1];
rz(1.6383621) q[1];
rz(-pi) q[2];
rz(2.2947512) q[3];
sx q[3];
rz(-1.2961642) q[3];
sx q[3];
rz(-1.3534897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.148823) q[2];
sx q[2];
rz(-1.6042446) q[2];
sx q[2];
rz(2.0796622) q[2];
rz(-0.81260931) q[3];
sx q[3];
rz(-1.142623) q[3];
sx q[3];
rz(-2.699471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84520311) q[0];
sx q[0];
rz(-2.3136105) q[0];
sx q[0];
rz(-1.4985341) q[0];
rz(2.8614112) q[1];
sx q[1];
rz(-1.2138841) q[1];
sx q[1];
rz(-0.79634204) q[1];
rz(-1.5308357) q[2];
sx q[2];
rz(-1.9794977) q[2];
sx q[2];
rz(-0.13499311) q[2];
rz(0.82874684) q[3];
sx q[3];
rz(-1.8846877) q[3];
sx q[3];
rz(1.7170513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
