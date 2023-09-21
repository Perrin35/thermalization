OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(2.467632) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622409) q[0];
sx q[0];
rz(-1.9277713) q[0];
sx q[0];
rz(0.39874052) q[0];
rz(-pi) q[1];
rz(1.7180213) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(2.0369612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4218581) q[1];
sx q[1];
rz(-2.0761479) q[1];
sx q[1];
rz(1.964633) q[1];
rz(-pi) q[2];
rz(1.854033) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42694416) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(-0.0077008458) q[0];
x q[1];
rz(-3.0395095) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(-1.724147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97034772) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(2.99511) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0262358) q[3];
sx q[3];
rz(-0.028209837) q[3];
sx q[3];
rz(-1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(-2.9648932) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(-1.7938991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189258) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(2.6618883) q[0];
x q[1];
rz(1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(-2.3929838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8828391) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-3.019481) q[1];
rz(-pi) q[2];
rz(-1.3860116) q[3];
sx q[3];
rz(-1.6596969) q[3];
sx q[3];
rz(0.88528663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72258654) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(-1.0976085) q[0];
x q[1];
rz(-0.16673659) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(-2.2819448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5986961) q[1];
sx q[1];
rz(-2.9196432) q[1];
sx q[1];
rz(1.0735682) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-0.6696223) q[3];
sx q[3];
rz(1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(0.66666493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71320888) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(1.3165738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(1.0416043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1866245) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(2.7093922) q[1];
x q[2];
rz(-2.2090962) q[3];
sx q[3];
rz(-0.36358788) q[3];
sx q[3];
rz(0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(-1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996457) q[0];
sx q[0];
rz(-2.7149704) q[0];
sx q[0];
rz(-2.9645779) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7026664) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(2.211328) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.02009) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(2.9083088) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7641894) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(-1.5148439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-3.1121758) q[2];
rz(1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32719192) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(2.8776684) q[0];
rz(-pi) q[1];
rz(2.3235544) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4335945) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(-0.57452332) q[1];
rz(-pi) q[2];
rz(2.0931582) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(-0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-0.62057173) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(2.506315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3978183) q[0];
sx q[0];
rz(-1.9282856) q[0];
sx q[0];
rz(-2.9767569) q[0];
rz(-pi) q[1];
rz(-1.073457) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(1.2129601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5105671) q[1];
sx q[1];
rz(-0.67093611) q[1];
sx q[1];
rz(-2.7061694) q[1];
rz(0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(-2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-0.90973967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(-0.33466848) q[0];
rz(-0.34809525) q[2];
sx q[2];
rz(-1.3653792) q[2];
sx q[2];
rz(1.5945438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(2.9305305) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89197253) q[3];
sx q[3];
rz(-2.744305) q[3];
sx q[3];
rz(0.15115034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(-0.48661423) q[2];
rz(-3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(0.89458481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54297011) q[0];
sx q[0];
rz(-2.2036836) q[0];
sx q[0];
rz(-3.0892239) q[0];
rz(-pi) q[1];
rz(1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(1.4354401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2243832) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(1.9320095) q[1];
rz(-1.7747202) q[3];
sx q[3];
rz(-2.4813934) q[3];
sx q[3];
rz(0.064299671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(0.62906229) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(2.2330877) q[2];
sx q[2];
rz(-2.0178595) q[2];
sx q[2];
rz(3.0541228) q[2];
rz(-1.4150423) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];