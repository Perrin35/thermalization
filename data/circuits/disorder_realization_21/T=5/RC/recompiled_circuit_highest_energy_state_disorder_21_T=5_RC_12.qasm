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
rz(1.1462829) q[0];
sx q[0];
rz(-3.1132071) q[0];
sx q[0];
rz(2.1838768) q[0];
rz(0.43342844) q[1];
sx q[1];
rz(-1.537701) q[1];
sx q[1];
rz(2.5133207) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7341029) q[0];
sx q[0];
rz(-0.39034319) q[0];
sx q[0];
rz(2.8170308) q[0];
rz(2.2987709) q[2];
sx q[2];
rz(-2.5113238) q[2];
sx q[2];
rz(-2.3284019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1294259) q[1];
sx q[1];
rz(-2.2291227) q[1];
sx q[1];
rz(-0.40947394) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7946934) q[3];
sx q[3];
rz(-0.42169398) q[3];
sx q[3];
rz(1.2661915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6750703) q[2];
sx q[2];
rz(-2.1212981) q[2];
sx q[2];
rz(-0.38727078) q[2];
rz(-1.9668503) q[3];
sx q[3];
rz(-1.4459123) q[3];
sx q[3];
rz(0.89850473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77884498) q[0];
sx q[0];
rz(-1.400482) q[0];
sx q[0];
rz(1.269302) q[0];
rz(1.4461888) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(0.92099014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5808153) q[0];
sx q[0];
rz(-1.8660424) q[0];
sx q[0];
rz(-0.92714374) q[0];
rz(-pi) q[1];
rz(-0.66514449) q[2];
sx q[2];
rz(-1.7153622) q[2];
sx q[2];
rz(1.4276366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2650406) q[1];
sx q[1];
rz(-1.3807793) q[1];
sx q[1];
rz(-1.9273619) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1109741) q[3];
sx q[3];
rz(-0.22600442) q[3];
sx q[3];
rz(-1.6219559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15128073) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(-0.82378236) q[2];
rz(-1.524205) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(1.9302146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23564944) q[0];
sx q[0];
rz(-0.88053572) q[0];
sx q[0];
rz(-0.5916416) q[0];
rz(-0.94332424) q[1];
sx q[1];
rz(-1.0572409) q[1];
sx q[1];
rz(1.0234157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9236301) q[0];
sx q[0];
rz(-1.4537088) q[0];
sx q[0];
rz(-0.10770672) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.01565) q[2];
sx q[2];
rz(-1.890939) q[2];
sx q[2];
rz(-0.61540907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4893787) q[1];
sx q[1];
rz(-1.4783597) q[1];
sx q[1];
rz(0.031186179) q[1];
rz(0.68444201) q[3];
sx q[3];
rz(-1.8307643) q[3];
sx q[3];
rz(2.3459159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25980514) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(-1.3938495) q[2];
rz(1.4272089) q[3];
sx q[3];
rz(-0.87427679) q[3];
sx q[3];
rz(-0.29903308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77089906) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(2.4942177) q[0];
rz(-0.24230832) q[1];
sx q[1];
rz(-1.4360042) q[1];
sx q[1];
rz(1.4586331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1636617) q[0];
sx q[0];
rz(-2.7898698) q[0];
sx q[0];
rz(-0.39885421) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3056612) q[2];
sx q[2];
rz(-0.80785895) q[2];
sx q[2];
rz(0.69617803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9414166) q[1];
sx q[1];
rz(-1.5343018) q[1];
sx q[1];
rz(2.0564276) q[1];
rz(0.69497739) q[3];
sx q[3];
rz(-0.91013346) q[3];
sx q[3];
rz(-1.7373178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4316537) q[2];
sx q[2];
rz(-0.84448758) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(0.25105181) q[3];
sx q[3];
rz(-0.71728388) q[3];
sx q[3];
rz(2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5421903) q[0];
sx q[0];
rz(-2.0090065) q[0];
sx q[0];
rz(-2.9436924) q[0];
rz(1.2359765) q[1];
sx q[1];
rz(-2.7695152) q[1];
sx q[1];
rz(3.1370251) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2777956) q[0];
sx q[0];
rz(-2.5044873) q[0];
sx q[0];
rz(-0.4235913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0367583) q[2];
sx q[2];
rz(-1.4994086) q[2];
sx q[2];
rz(-2.57881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63268836) q[1];
sx q[1];
rz(-1.4211854) q[1];
sx q[1];
rz(2.8652111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6130535) q[3];
sx q[3];
rz(-2.740707) q[3];
sx q[3];
rz(-1.8141845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68134754) q[2];
sx q[2];
rz(-2.6474417) q[2];
sx q[2];
rz(-2.8013308) q[2];
rz(1.5334689) q[3];
sx q[3];
rz(-1.3932649) q[3];
sx q[3];
rz(-1.5197915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.915864) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(-1.1874143) q[0];
rz(1.4215218) q[1];
sx q[1];
rz(-2.7233796) q[1];
sx q[1];
rz(0.13030599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73459816) q[0];
sx q[0];
rz(-0.13665527) q[0];
sx q[0];
rz(1.1278485) q[0];
rz(-pi) q[1];
rz(-2.2537986) q[2];
sx q[2];
rz(-1.0947795) q[2];
sx q[2];
rz(2.7107875) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7728945) q[1];
sx q[1];
rz(-2.9737568) q[1];
sx q[1];
rz(3.1375727) q[1];
rz(3.1316787) q[3];
sx q[3];
rz(-1.8335206) q[3];
sx q[3];
rz(2.3309121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80798951) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(-1.0398593) q[2];
rz(-0.17397675) q[3];
sx q[3];
rz(-1.1287374) q[3];
sx q[3];
rz(0.66966331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60798764) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(-2.3002891) q[0];
rz(-1.816642) q[1];
sx q[1];
rz(-0.94580301) q[1];
sx q[1];
rz(-1.4680877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7588336) q[0];
sx q[0];
rz(-1.2099203) q[0];
sx q[0];
rz(1.1961236) q[0];
rz(-pi) q[1];
rz(3.1300342) q[2];
sx q[2];
rz(-2.0894139) q[2];
sx q[2];
rz(0.29588884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12097479) q[1];
sx q[1];
rz(-2.1281895) q[1];
sx q[1];
rz(1.0110028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2135336) q[3];
sx q[3];
rz(-1.8566086) q[3];
sx q[3];
rz(0.88102007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86772052) q[2];
sx q[2];
rz(-0.90110675) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(-1.0817179) q[3];
sx q[3];
rz(-0.41926256) q[3];
sx q[3];
rz(0.7209512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69407392) q[0];
sx q[0];
rz(-0.17386359) q[0];
sx q[0];
rz(-0.41282594) q[0];
rz(1.5326356) q[1];
sx q[1];
rz(-2.2282579) q[1];
sx q[1];
rz(-2.434381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21239534) q[0];
sx q[0];
rz(-1.2629422) q[0];
sx q[0];
rz(-2.3677127) q[0];
rz(-pi) q[1];
rz(2.9712255) q[2];
sx q[2];
rz(-2.3922386) q[2];
sx q[2];
rz(0.10766497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2688931) q[1];
sx q[1];
rz(-1.2352691) q[1];
sx q[1];
rz(2.5771532) q[1];
rz(1.0319581) q[3];
sx q[3];
rz(-1.4488359) q[3];
sx q[3];
rz(-1.7820047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1327208) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(-0.51187619) q[2];
rz(-0.68073186) q[3];
sx q[3];
rz(-1.3635037) q[3];
sx q[3];
rz(2.9752922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7452288) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(0.45528278) q[0];
rz(-1.0633172) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(-0.071970073) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44350478) q[0];
sx q[0];
rz(-0.035497276) q[0];
sx q[0];
rz(-2.7867691) q[0];
rz(-pi) q[1];
rz(-2.7006702) q[2];
sx q[2];
rz(-1.5917695) q[2];
sx q[2];
rz(1.4628589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0177415) q[1];
sx q[1];
rz(-1.2203274) q[1];
sx q[1];
rz(0.88092901) q[1];
rz(-pi) q[2];
x q[2];
rz(1.23986) q[3];
sx q[3];
rz(-1.8171105) q[3];
sx q[3];
rz(-0.34231753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7887743) q[2];
sx q[2];
rz(-2.7249551) q[2];
sx q[2];
rz(-0.62208661) q[2];
rz(2.3173053) q[3];
sx q[3];
rz(-1.0930748) q[3];
sx q[3];
rz(2.2881959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033009919) q[0];
sx q[0];
rz(-2.0350631) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(-1.7225601) q[1];
sx q[1];
rz(-2.6585237) q[1];
sx q[1];
rz(2.5249265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4224859) q[0];
sx q[0];
rz(-2.2107443) q[0];
sx q[0];
rz(-1.5404705) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5316467) q[2];
sx q[2];
rz(-1.8862572) q[2];
sx q[2];
rz(-1.1300398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2076555) q[1];
sx q[1];
rz(-0.21156921) q[1];
sx q[1];
rz(0.013222828) q[1];
rz(-pi) q[2];
rz(-1.0000049) q[3];
sx q[3];
rz(-0.94729187) q[3];
sx q[3];
rz(-2.4542261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0509433) q[2];
sx q[2];
rz(-1.1530777) q[2];
sx q[2];
rz(1.0672182) q[2];
rz(2.0459335) q[3];
sx q[3];
rz(-1.9039543) q[3];
sx q[3];
rz(0.9637951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6296366) q[0];
sx q[0];
rz(-2.1558599) q[0];
sx q[0];
rz(2.070367) q[0];
rz(1.4818954) q[1];
sx q[1];
rz(-1.0922468) q[1];
sx q[1];
rz(-2.560871) q[1];
rz(2.8790375) q[2];
sx q[2];
rz(-1.5468183) q[2];
sx q[2];
rz(0.47894947) q[2];
rz(0.4023424) q[3];
sx q[3];
rz(-0.33807031) q[3];
sx q[3];
rz(-2.8344179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
