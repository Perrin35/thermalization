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
rz(-1.988451) q[0];
sx q[0];
rz(-2.3260131) q[0];
sx q[0];
rz(-2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(4.9795436) q[1];
sx q[1];
rz(10.995168) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2822097) q[0];
sx q[0];
rz(-2.4042685) q[0];
sx q[0];
rz(2.5928549) q[0];
rz(-0.067367359) q[2];
sx q[2];
rz(-1.172003) q[2];
sx q[2];
rz(0.68902868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0834404) q[1];
sx q[1];
rz(-1.6279164) q[1];
sx q[1];
rz(-1.2099464) q[1];
rz(0.060293555) q[3];
sx q[3];
rz(-1.2620842) q[3];
sx q[3];
rz(-1.2293881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5300753) q[2];
sx q[2];
rz(-0.0081743058) q[2];
sx q[2];
rz(0.44069904) q[2];
rz(0.079744451) q[3];
sx q[3];
rz(-3.1414882) q[3];
sx q[3];
rz(1.124148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663064) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(-2.9735907) q[0];
rz(-0.02027823) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(1.5365938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0609866) q[0];
sx q[0];
rz(-1.8053375) q[0];
sx q[0];
rz(-1.7936262) q[0];
x q[1];
rz(-1.3290908) q[2];
sx q[2];
rz(-0.010420825) q[2];
sx q[2];
rz(-1.3543606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3520842) q[1];
sx q[1];
rz(-1.5696973) q[1];
sx q[1];
rz(1.5754682) q[1];
rz(-0.042293799) q[3];
sx q[3];
rz(-1.5618069) q[3];
sx q[3];
rz(0.90153722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-2.2296843) q[2];
sx q[2];
rz(1.7339285) q[2];
rz(-1.0450854) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(0.27200562) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8318091) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(-0.56104863) q[0];
rz(-2.8647515) q[1];
sx q[1];
rz(-0.012877348) q[1];
sx q[1];
rz(1.8337839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9914068) q[0];
sx q[0];
rz(-1.5298109) q[0];
sx q[0];
rz(-1.2738373) q[0];
rz(3.1342034) q[2];
sx q[2];
rz(-1.5708013) q[2];
sx q[2];
rz(-2.6342692) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11266358) q[1];
sx q[1];
rz(-1.512292) q[1];
sx q[1];
rz(-0.99716352) q[1];
rz(-pi) q[2];
rz(0.8933634) q[3];
sx q[3];
rz(-1.2033389) q[3];
sx q[3];
rz(1.2585672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7967367) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(2.5522088) q[2];
rz(-1.899259) q[3];
sx q[3];
rz(-3.1292249) q[3];
sx q[3];
rz(-1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068483343) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(0.033500813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.861167) q[0];
sx q[0];
rz(-0.67641947) q[0];
sx q[0];
rz(-1.2926213) q[0];
rz(-0.11963614) q[2];
sx q[2];
rz(-1.5663317) q[2];
sx q[2];
rz(1.9890832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9814947) q[1];
sx q[1];
rz(-1.3019053) q[1];
sx q[1];
rz(3.1310496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93422814) q[3];
sx q[3];
rz(-2.9381972) q[3];
sx q[3];
rz(2.0980841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.032430705) q[2];
sx q[2];
rz(-3.1350632) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.2700891) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700767) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(0.44684967) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(-1.7284547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2010445) q[0];
sx q[0];
rz(-1.5371635) q[0];
sx q[0];
rz(1.3399009) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1485708) q[2];
sx q[2];
rz(-1.7686426) q[2];
sx q[2];
rz(-2.5289218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.079559673) q[1];
sx q[1];
rz(-1.6158982) q[1];
sx q[1];
rz(-1.5212039) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2893612) q[3];
sx q[3];
rz(-1.7858943) q[3];
sx q[3];
rz(0.18751442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3073005) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(2.6429122) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(-2.5448866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805098) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(2.7951796) q[0];
rz(0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(2.3880889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2997871) q[0];
sx q[0];
rz(-2.8827169) q[0];
sx q[0];
rz(2.4095834) q[0];
x q[1];
rz(0.97917231) q[2];
sx q[2];
rz(-3.0074928) q[2];
sx q[2];
rz(-0.025452415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1767039) q[1];
sx q[1];
rz(-2.2060925) q[1];
sx q[1];
rz(0.52365644) q[1];
rz(1.7338792) q[3];
sx q[3];
rz(-1.4797416) q[3];
sx q[3];
rz(-2.1084821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5715013) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(1.6174779) q[2];
rz(-0.12024719) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(1.6809173) q[1];
sx q[1];
rz(-2.2082081) q[1];
sx q[1];
rz(-0.077300765) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6344096) q[0];
sx q[0];
rz(-1.612525) q[0];
sx q[0];
rz(-0.001464837) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5796698) q[2];
sx q[2];
rz(-1.5755782) q[2];
sx q[2];
rz(1.3591131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12880023) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(2.7526593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3568001) q[3];
sx q[3];
rz(-1.4460751) q[3];
sx q[3];
rz(2.7217442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7920502) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(2.1816317) q[2];
rz(0.32944426) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9462117) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(-0.10928133) q[0];
rz(0.37336135) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(1.9111309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41540158) q[0];
sx q[0];
rz(-2.2080732) q[0];
sx q[0];
rz(-2.3660052) q[0];
x q[1];
rz(0.76704278) q[2];
sx q[2];
rz(-2.8708176) q[2];
sx q[2];
rz(-2.3363638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8459863) q[1];
sx q[1];
rz(-1.5822268) q[1];
sx q[1];
rz(1.5039526) q[1];
rz(-0.30970311) q[3];
sx q[3];
rz(-0.94325698) q[3];
sx q[3];
rz(2.5554339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5750778) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(1.396842) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0777271) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(2.5685837) q[0];
rz(-0.30814463) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(-1.0073957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1290179) q[0];
sx q[0];
rz(-1.5652302) q[0];
sx q[0];
rz(1.6773423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2635055) q[2];
sx q[2];
rz(-2.9973534) q[2];
sx q[2];
rz(-2.7968614) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17698174) q[1];
sx q[1];
rz(-1.4876502) q[1];
sx q[1];
rz(-0.040772922) q[1];
rz(-pi) q[2];
x q[2];
rz(0.013434826) q[3];
sx q[3];
rz(-1.6009496) q[3];
sx q[3];
rz(2.9363901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8241626) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(2.7516348) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(2.2700229) q[1];
sx q[1];
rz(-1.8337245) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594306) q[0];
sx q[0];
rz(-2.5100187) q[0];
sx q[0];
rz(2.6078014) q[0];
rz(-pi) q[1];
rz(1.5266085) q[2];
sx q[2];
rz(-0.95643294) q[2];
sx q[2];
rz(-0.0690661) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20032665) q[1];
sx q[1];
rz(-1.8706053) q[1];
sx q[1];
rz(2.7977562) q[1];
rz(-pi) q[2];
rz(-1.6127729) q[3];
sx q[3];
rz(-1.4577565) q[3];
sx q[3];
rz(2.6338995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5747052) q[2];
sx q[2];
rz(-3.0988099) q[2];
sx q[2];
rz(-0.031878397) q[2];
rz(2.3632862) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(-2.8460898) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185709) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(0.12693916) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(1.610582) q[2];
sx q[2];
rz(-3.0018158) q[2];
sx q[2];
rz(-2.9374585) q[2];
rz(-0.89613468) q[3];
sx q[3];
rz(-1.5265161) q[3];
sx q[3];
rz(2.0731887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
