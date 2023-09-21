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
rz(-0.67396069) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062382467) q[0];
sx q[0];
rz(-1.198472) q[0];
sx q[0];
rz(-1.9553493) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27172471) q[2];
sx q[2];
rz(-0.50479111) q[2];
sx q[2];
rz(1.4128078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4218581) q[1];
sx q[1];
rz(-2.0761479) q[1];
sx q[1];
rz(-1.964633) q[1];
rz(-pi) q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1439046) q[0];
sx q[0];
rz(-1.5630957) q[0];
sx q[0];
rz(1.5776004) q[0];
rz(-pi) q[1];
rz(-3.0395095) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(1.4174457) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56541601) q[1];
sx q[1];
rz(-1.4285354) q[1];
sx q[1];
rz(-1.3284645) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1291817) q[3];
sx q[3];
rz(-1.5454626) q[3];
sx q[3];
rz(-2.2001571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-0.97989782) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(-1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-2.9648932) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(-1.7938991) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603148) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(-0.48849948) q[0];
rz(3.0323896) q[2];
sx q[2];
rz(-2.9401527) q[2];
sx q[2];
rz(0.86004721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(0.12211166) q[1];
x q[2];
rz(1.3860116) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(-2.256306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692114) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4190061) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(-2.0439842) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16673659) q[2];
sx q[2];
rz(-2.1519289) q[2];
sx q[2];
rz(2.2819448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1064062) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(3.0343642) q[1];
rz(-pi) q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-0.6696223) q[3];
sx q[3];
rz(1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(0.60194683) q[2];
rz(-2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
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
rz(-2.4749277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71320888) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(1.3165738) q[0];
rz(2.1140852) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(1.0416043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40068377) q[1];
sx q[1];
rz(-2.2212257) q[1];
sx q[1];
rz(1.9294192) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9186451) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(-1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.5302352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2419469) q[0];
sx q[0];
rz(-2.7149704) q[0];
sx q[0];
rz(-0.17701478) q[0];
rz(-2.5189581) q[2];
sx q[2];
rz(-2.6175675) q[2];
sx q[2];
rz(0.084409075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64165243) q[1];
sx q[1];
rz(-1.3556726) q[1];
sx q[1];
rz(-1.9745449) q[1];
rz(-pi) q[2];
rz(-0.050614428) q[3];
sx q[3];
rz(-1.8236056) q[3];
sx q[3];
rz(1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683838) q[0];
sx q[0];
rz(-1.3084992) q[0];
sx q[0];
rz(1.6845076) q[0];
x q[1];
rz(-2.4099318) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(2.4351956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43607298) q[1];
sx q[1];
rz(-1.1735859) q[1];
sx q[1];
rz(-0.7049837) q[1];
rz(-pi) q[2];
rz(2.8227311) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(1.5368411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(-2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1875293) q[0];
sx q[0];
rz(-0.39217338) q[0];
sx q[0];
rz(-1.9847045) q[0];
rz(0.17417553) q[2];
sx q[2];
rz(-1.0798228) q[2];
sx q[2];
rz(2.8665286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(2.5177588) q[1];
rz(-pi) q[2];
rz(-0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(-3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(3.1414202) q[2];
rz(-2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-0.85787684) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37659392) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(-1.8565208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-2.6532432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(-0.21106212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.886456) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(2.6549784) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.2238812) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(2.2470078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986225) q[0];
sx q[0];
rz(-0.93790903) q[0];
sx q[0];
rz(0.052368725) q[0];
x q[1];
rz(-1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(2.6596136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67748653) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(-3.0739215) q[1];
x q[2];
rz(1.7747202) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(-3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(-1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(0.62906229) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(-3.0929052) q[3];
sx q[3];
rz(-1.2702474) q[3];
sx q[3];
rz(-2.4297759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
