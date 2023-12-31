OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(-1.8859099) q[0];
sx q[0];
rz(0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(-2.9405138) q[1];
sx q[1];
rz(-0.091436401) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(-0.63011516) q[0];
x q[1];
rz(2.6684127) q[2];
sx q[2];
rz(-2.8521529) q[2];
sx q[2];
rz(-0.48130408) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(1.9074349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9760194) q[3];
sx q[3];
rz(-1.752749) q[3];
sx q[3];
rz(0.44427696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66449195) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-2.9512761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7305671) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(-0.21582027) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9637738) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(2.9289392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39365921) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(-3.1118803) q[1];
x q[2];
rz(-2.4137647) q[3];
sx q[3];
rz(-1.9350633) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.2480199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(0.16288217) q[0];
rz(-pi) q[1];
rz(1.6349995) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.9931672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(-1.6559421) q[1];
rz(-0.30001254) q[3];
sx q[3];
rz(-0.94167751) q[3];
sx q[3];
rz(2.574207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(1.8444555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545249) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(-1.7012419) q[0];
x q[1];
rz(-2.9739431) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(2.167785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0399196) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(-2.4458829) q[1];
x q[2];
rz(-1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(1.0877346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.8937768) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
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
rz(-2.8844937) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(2.2280072) q[0];
rz(2.2141453) q[2];
sx q[2];
rz(-0.38380917) q[2];
sx q[2];
rz(-0.62188934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51902276) q[1];
sx q[1];
rz(-1.6791324) q[1];
sx q[1];
rz(0.83801724) q[1];
x q[2];
rz(1.250508) q[3];
sx q[3];
rz(-2.3158145) q[3];
sx q[3];
rz(-2.461493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(2.7048236) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540928) q[0];
sx q[0];
rz(-2.0797074) q[0];
sx q[0];
rz(-0.66977588) q[0];
rz(2.2366441) q[2];
sx q[2];
rz(-1.6709423) q[2];
sx q[2];
rz(0.03749321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.9532277) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(0.26785775) q[1];
rz(-2.3267641) q[3];
sx q[3];
rz(-2.1216672) q[3];
sx q[3];
rz(1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(2.4760822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367046) q[0];
sx q[0];
rz(-2.8993653) q[0];
sx q[0];
rz(-1.5664943) q[0];
rz(-2.9930816) q[2];
sx q[2];
rz(-0.87887895) q[2];
sx q[2];
rz(-1.6910545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.682714) q[1];
sx q[1];
rz(-0.51358089) q[1];
sx q[1];
rz(-2.0290124) q[1];
rz(-pi) q[2];
rz(2.9969278) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(-0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539464) q[0];
sx q[0];
rz(-2.7641312) q[0];
sx q[0];
rz(1.1520555) q[0];
rz(0.45192265) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-2.506633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6845477) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(-0.12462516) q[1];
x q[2];
rz(1.8975774) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(-2.9941878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.609628) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(2.4832446) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(0.13959612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0567719) q[0];
sx q[0];
rz(-0.064282566) q[0];
sx q[0];
rz(-1.2540741) q[0];
x q[1];
rz(0.058768674) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(1.2440484) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82397205) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(3.0082263) q[1];
rz(-pi) q[2];
rz(-3.0890907) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(-2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(1.2963294) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(-0.18558003) q[0];
rz(-1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(-1.0111965) q[0];
rz(-0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-0.05664209) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85877307) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(-1.8358843) q[1];
x q[2];
rz(0.0049403355) q[3];
sx q[3];
rz(-2.2757109) q[3];
sx q[3];
rz(2.4019394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(2.623693) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26375297) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-0.18763018) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-0.27992237) q[2];
sx q[2];
rz(-1.7114077) q[2];
sx q[2];
rz(-0.68243295) q[2];
rz(-2.256176) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
