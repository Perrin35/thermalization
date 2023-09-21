OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(2.8136301) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5782195) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(2.5114775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47317998) q[2];
sx q[2];
rz(-2.8521529) q[2];
sx q[2];
rz(0.48130408) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58127296) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(-1.2341577) q[1];
x q[2];
rz(0.21858469) q[3];
sx q[3];
rz(-0.9871452) q[3];
sx q[3];
rz(2.1368795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(2.8664355) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(2.9512761) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0424461) q[0];
sx q[0];
rz(-0.45396807) q[0];
sx q[0];
rz(1.1048261) q[0];
rz(-2.6066577) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(-1.5869706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1574993) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(-2.2929272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.042949) q[3];
sx q[3];
rz(-2.241579) q[3];
sx q[3];
rz(0.19876476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(2.9787105) q[0];
x q[1];
rz(-0.2072316) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(0.40916967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4856845) q[1];
sx q[1];
rz(-0.67042065) q[1];
sx q[1];
rz(-0.1078492) q[1];
x q[2];
rz(2.2218024) q[3];
sx q[3];
rz(-1.3295104) q[3];
sx q[3];
rz(1.1834708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.9281663) q[2];
rz(-0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(1.8444555) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521833) q[0];
sx q[0];
rz(-1.7011189) q[0];
sx q[0];
rz(-3.0979063) q[0];
rz(-2.0879891) q[2];
sx q[2];
rz(-0.32978168) q[2];
sx q[2];
rz(-1.5151378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2343826) q[1];
sx q[1];
rz(-1.0115336) q[1];
sx q[1];
rz(0.84749605) q[1];
rz(1.3912195) q[3];
sx q[3];
rz(-1.5571801) q[3];
sx q[3];
rz(2.7333958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(0.91919351) q[2];
rz(0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.325266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(-1.7153046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.250524) q[0];
sx q[0];
rz(-2.2254125) q[0];
sx q[0];
rz(-3.0380681) q[0];
x q[1];
rz(1.8833141) q[2];
sx q[2];
rz(-1.3442355) q[2];
sx q[2];
rz(0.3414008) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1868185) q[1];
sx q[1];
rz(-0.84328077) q[1];
sx q[1];
rz(0.14528841) q[1];
rz(-0.32894965) q[3];
sx q[3];
rz(-0.79862404) q[3];
sx q[3];
rz(-2.9165099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8067236) q[0];
sx q[0];
rz(-0.81672251) q[0];
sx q[0];
rz(2.4094894) q[0];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(-1.6598998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39396573) q[1];
sx q[1];
rz(-2.6193301) q[1];
sx q[1];
rz(-1.0737435) q[1];
x q[2];
rz(2.3267641) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(0.56419939) q[2];
rz(0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77008477) q[0];
sx q[0];
rz(-1.5718282) q[0];
sx q[0];
rz(-1.8130215) q[0];
rz(0.87343563) q[2];
sx q[2];
rz(-1.6849815) q[2];
sx q[2];
rz(0.21542491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.682714) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(2.0290124) q[1];
x q[2];
rz(-2.1036759) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(-2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-0.19259024) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(0.15983454) q[0];
x q[1];
rz(2.68967) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-0.6349596) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0847816) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(2.0463498) q[1];
rz(-pi) q[2];
rz(-0.87266463) q[3];
sx q[3];
rz(-0.41655311) q[3];
sx q[3];
rz(2.0779028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(1.2540741) q[0];
rz(-pi) q[1];
rz(2.326968) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(2.8551561) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(2.6194797) q[1];
rz(1.6406052) q[3];
sx q[3];
rz(-2.4959292) q[3];
sx q[3];
rz(2.0696236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6167986) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(0.087879114) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907222) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(-2.0275293) q[0];
x q[1];
rz(2.2514258) q[2];
sx q[2];
rz(-0.3496799) q[2];
sx q[2];
rz(-2.4868951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85877307) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(-1.3057083) q[1];
x q[2];
rz(-0.0049403355) q[3];
sx q[3];
rz(-2.2757109) q[3];
sx q[3];
rz(0.73965328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(0.55220848) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.26375297) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(0.4734584) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(0.88541661) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];