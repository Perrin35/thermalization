OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1333753) q[0];
sx q[0];
rz(-1.8741338) q[0];
sx q[0];
rz(-3.128669) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(2.0838783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81257129) q[0];
sx q[0];
rz(-2.2264535) q[0];
sx q[0];
rz(1.283487) q[0];
x q[1];
rz(-0.51527649) q[2];
sx q[2];
rz(-0.84723398) q[2];
sx q[2];
rz(1.1267337) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8243421) q[1];
sx q[1];
rz(-1.9320654) q[1];
sx q[1];
rz(0.74202219) q[1];
rz(-pi) q[2];
rz(0.042594508) q[3];
sx q[3];
rz(-1.9622246) q[3];
sx q[3];
rz(3.1393876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5554123) q[2];
sx q[2];
rz(-1.5988028) q[2];
sx q[2];
rz(-3.0482698) q[2];
rz(0.020545067) q[3];
sx q[3];
rz(-13/(16*pi)) q[3];
sx q[3];
rz(-1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071863197) q[0];
sx q[0];
rz(-1.3988031) q[0];
sx q[0];
rz(2.3139957) q[0];
rz(0.96356511) q[1];
sx q[1];
rz(-1.6160485) q[1];
sx q[1];
rz(2.4049984) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73294357) q[0];
sx q[0];
rz(-1.2970862) q[0];
sx q[0];
rz(0.078775551) q[0];
rz(-pi) q[1];
rz(2.5649293) q[2];
sx q[2];
rz(-1.0500337) q[2];
sx q[2];
rz(1.3530089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40562661) q[1];
sx q[1];
rz(-2.7752004) q[1];
sx q[1];
rz(-0.97586164) q[1];
rz(-pi) q[2];
rz(1.5932036) q[3];
sx q[3];
rz(-1.4364035) q[3];
sx q[3];
rz(-2.973864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57614148) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(0.60892504) q[3];
sx q[3];
rz(-2.3600793) q[3];
sx q[3];
rz(1.5003381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(2.3420077) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(2.0094357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6511667) q[0];
sx q[0];
rz(-1.4803783) q[0];
sx q[0];
rz(-1.7866485) q[0];
x q[1];
rz(-3.0614616) q[2];
sx q[2];
rz(-1.5485768) q[2];
sx q[2];
rz(-0.20168389) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70402789) q[1];
sx q[1];
rz(-1.6482185) q[1];
sx q[1];
rz(0.67005007) q[1];
rz(-pi) q[2];
rz(1.5270832) q[3];
sx q[3];
rz(-0.96230405) q[3];
sx q[3];
rz(0.19245806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-2.5567882) q[2];
sx q[2];
rz(-1.7542138) q[2];
rz(-0.34902188) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(-0.52946985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.741852) q[0];
sx q[0];
rz(-0.96804237) q[0];
sx q[0];
rz(-0.35811785) q[0];
rz(-2.070836) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(-1.4395641) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5191906) q[0];
sx q[0];
rz(-0.82525142) q[0];
sx q[0];
rz(-0.36143266) q[0];
rz(-pi) q[1];
rz(0.49825495) q[2];
sx q[2];
rz(-0.94211218) q[2];
sx q[2];
rz(1.2815338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6144952) q[1];
sx q[1];
rz(-2.1367837) q[1];
sx q[1];
rz(-2.4742592) q[1];
rz(-pi) q[2];
rz(2.8791588) q[3];
sx q[3];
rz(-2.35733) q[3];
sx q[3];
rz(-1.5465496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4293356) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(-2.5168391) q[2];
rz(1.6338927) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(-2.0190575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37609142) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(-1.435185) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(2.3211839) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89264458) q[0];
sx q[0];
rz(-1.2995509) q[0];
sx q[0];
rz(0.091288996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7960288) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(2.1175543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0972406) q[1];
sx q[1];
rz(-0.28438452) q[1];
sx q[1];
rz(-0.76459717) q[1];
rz(-2.5518604) q[3];
sx q[3];
rz(-1.5431649) q[3];
sx q[3];
rz(-3.0363135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72941214) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(-2.6030276) q[2];
rz(1.4696848) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84022123) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(0.090601623) q[0];
rz(-2.7719356) q[1];
sx q[1];
rz(-1.5934817) q[1];
sx q[1];
rz(2.7634117) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2205296) q[0];
sx q[0];
rz(-2.605633) q[0];
sx q[0];
rz(-0.45951636) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43439602) q[2];
sx q[2];
rz(-3.042331) q[2];
sx q[2];
rz(2.1289265) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82434618) q[1];
sx q[1];
rz(-2.6326944) q[1];
sx q[1];
rz(-0.25516487) q[1];
x q[2];
rz(1.7615116) q[3];
sx q[3];
rz(-1.4927974) q[3];
sx q[3];
rz(1.5275148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52046627) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(1.8488047) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(1.7267936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3066278) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(0.90245885) q[0];
rz(0.2392256) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(-0.36144027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1162565) q[0];
sx q[0];
rz(-0.5098719) q[0];
sx q[0];
rz(-1.8007397) q[0];
rz(-pi) q[1];
rz(-3.0143033) q[2];
sx q[2];
rz(-1.5152928) q[2];
sx q[2];
rz(-1.0148359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8351961) q[1];
sx q[1];
rz(-1.3521242) q[1];
sx q[1];
rz(1.1347215) q[1];
x q[2];
rz(-0.33935762) q[3];
sx q[3];
rz(-1.7202783) q[3];
sx q[3];
rz(-0.71163346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0685737) q[2];
sx q[2];
rz(-1.8014427) q[2];
sx q[2];
rz(-2.2646591) q[2];
rz(0.98313156) q[3];
sx q[3];
rz(-1.5777595) q[3];
sx q[3];
rz(0.94299281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.8125732) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(0.83874291) q[0];
rz(-2.4328531) q[1];
sx q[1];
rz(-0.63116169) q[1];
sx q[1];
rz(2.2355524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748286) q[0];
sx q[0];
rz(-1.5918918) q[0];
sx q[0];
rz(2.4167908) q[0];
rz(-2.110741) q[2];
sx q[2];
rz(-0.40274061) q[2];
sx q[2];
rz(2.7583721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.816765) q[1];
sx q[1];
rz(-0.93597368) q[1];
sx q[1];
rz(-0.88714182) q[1];
rz(2.1791547) q[3];
sx q[3];
rz(-1.7791074) q[3];
sx q[3];
rz(2.359451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9771542) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(-1.1158811) q[2];
rz(0.62853938) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(-0.61454296) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81166613) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(-0.16127583) q[0];
rz(2.2992086) q[1];
sx q[1];
rz(-1.7120275) q[1];
sx q[1];
rz(-2.9715723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1415213) q[0];
sx q[0];
rz(-1.562644) q[0];
sx q[0];
rz(3.1233112) q[0];
rz(-2.7964994) q[2];
sx q[2];
rz(-2.6775914) q[2];
sx q[2];
rz(-2.4008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57178088) q[1];
sx q[1];
rz(-1.7649635) q[1];
sx q[1];
rz(-1.1295736) q[1];
rz(-pi) q[2];
rz(-0.028548553) q[3];
sx q[3];
rz(-1.1880837) q[3];
sx q[3];
rz(0.0023502758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1539479) q[2];
sx q[2];
rz(-1.597007) q[2];
sx q[2];
rz(1.4863996) q[2];
rz(0.060128309) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.7000343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104367) q[0];
sx q[0];
rz(-0.031351723) q[0];
sx q[0];
rz(-1.4309058) q[0];
rz(0.36262861) q[1];
sx q[1];
rz(-1.5321833) q[1];
sx q[1];
rz(1.8261725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.423374) q[0];
sx q[0];
rz(-0.43180433) q[0];
sx q[0];
rz(-0.81599094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9489297) q[2];
sx q[2];
rz(-1.6844498) q[2];
sx q[2];
rz(0.092157539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23116701) q[1];
sx q[1];
rz(-0.98955123) q[1];
sx q[1];
rz(1.890593) q[1];
x q[2];
rz(0.4591367) q[3];
sx q[3];
rz(-0.92184421) q[3];
sx q[3];
rz(0.54455633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0843087) q[2];
sx q[2];
rz(-1.5466651) q[2];
sx q[2];
rz(2.9810737) q[2];
rz(-2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01874825) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(1.979076) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(3.054259) q[2];
sx q[2];
rz(-1.4764016) q[2];
sx q[2];
rz(1.7511677) q[2];
rz(2.3707262) q[3];
sx q[3];
rz(-1.1419747) q[3];
sx q[3];
rz(2.9621073) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
