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
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51973242) q[0];
sx q[0];
rz(-0.73017263) q[0];
sx q[0];
rz(-0.61868389) q[0];
rz(2.882471) q[2];
sx q[2];
rz(-1.7012351) q[2];
sx q[2];
rz(-2.5082617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2135703) q[1];
sx q[1];
rz(-1.3334647) q[1];
sx q[1];
rz(0.084753239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8880635) q[3];
sx q[3];
rz(-2.5228365) q[3];
sx q[3];
rz(1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-2.9512761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961401) q[0];
sx q[0];
rz(-1.7691233) q[0];
sx q[0];
rz(-1.9818927) q[0];
rz(2.6066577) q[2];
sx q[2];
rz(-1.9133647) q[2];
sx q[2];
rz(-1.5869706) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9840934) q[1];
sx q[1];
rz(-1.5930953) q[1];
sx q[1];
rz(0.84866546) q[1];
rz(-pi) q[2];
x q[2];
rz(2.042949) q[3];
sx q[3];
rz(-0.90001366) q[3];
sx q[3];
rz(0.19876476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(-1.2480199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67278636) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(2.9787105) q[0];
rz(0.2072316) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(-0.40916967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51860147) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(-1.6559421) q[1];
rz(-pi) q[2];
rz(1.1850584) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(0.08337534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(-3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064917795) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(1.8924367) q[0];
rz(-pi) q[1];
x q[1];
rz(1.859971) q[2];
sx q[2];
rz(-1.4099858) q[2];
sx q[2];
rz(-0.5493872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(-0.69570978) q[1];
x q[2];
rz(1.6468871) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(-2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.7153046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25709897) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(2.2280072) q[0];
rz(-0.9274474) q[2];
sx q[2];
rz(-2.7577835) q[2];
sx q[2];
rz(0.62188934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9702455) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(-1.7319748) q[1];
rz(0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(-2.9165099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(-0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(2.4018535) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(0.57156634) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7874998) q[0];
sx q[0];
rz(-1.0618853) q[0];
sx q[0];
rz(-0.66977588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4095441) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(-1.4816928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.188365) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(2.8737349) q[1];
x q[2];
rz(-0.84047517) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(-0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(0.71227658) q[0];
rz(-0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8004566) q[0];
sx q[0];
rz(-1.3285713) q[0];
sx q[0];
rz(-0.0010629396) q[0];
rz(-pi) q[1];
rz(2.9930816) q[2];
sx q[2];
rz(-0.87887895) q[2];
sx q[2];
rz(-1.4505381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1674541) q[1];
sx q[1];
rz(-2.0270837) q[1];
sx q[1];
rz(-0.2445226) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0379167) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(-0.52475196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(2.802882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324547) q[0];
sx q[0];
rz(-1.7212241) q[0];
sx q[0];
rz(1.223279) q[0];
x q[1];
rz(1.8024826) q[2];
sx q[2];
rz(-1.1294239) q[2];
sx q[2];
rz(1.0362792) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72570669) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.8094256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2440153) q[3];
sx q[3];
rz(-1.3076926) q[3];
sx q[3];
rz(0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(0.65834808) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-0.13959612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(-1.8875185) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5084969) q[2];
sx q[2];
rz(-2.3844516) q[2];
sx q[2];
rz(-1.8119259) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3176206) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(-0.13336639) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052501909) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(-2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(0.33111462) q[2];
rz(-2.3838499) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.4846444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7916242) q[0];
sx q[0];
rz(-2.0098643) q[0];
sx q[0];
rz(-0.29859782) q[0];
x q[1];
rz(-1.8469641) q[2];
sx q[2];
rz(-1.7880926) q[2];
sx q[2];
rz(1.5751788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(1.9468716) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86587571) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-2.6681343) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(2.3579303) q[3];
sx q[3];
rz(-0.85859921) q[3];
sx q[3];
rz(2.8006299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
