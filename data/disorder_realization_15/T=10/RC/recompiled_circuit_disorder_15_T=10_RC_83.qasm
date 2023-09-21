OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29785922) q[0];
sx q[0];
rz(-2.5279186) q[0];
sx q[0];
rz(2.4224129) q[0];
rz(-1.7742046) q[1];
sx q[1];
rz(-2.8957638) q[1];
sx q[1];
rz(-2.153102) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0713232) q[0];
sx q[0];
rz(-1.9914314) q[0];
sx q[0];
rz(-0.063739901) q[0];
rz(-0.97889401) q[2];
sx q[2];
rz(-1.7040633) q[2];
sx q[2];
rz(-2.8507581) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5060252) q[1];
sx q[1];
rz(-2.4534181) q[1];
sx q[1];
rz(0.78189214) q[1];
rz(-pi) q[2];
rz(-3.1036166) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(-1.2167041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6926379) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.9521936) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-2.7785595) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(1.2526858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38329268) q[0];
sx q[0];
rz(-2.9581666) q[0];
sx q[0];
rz(2.6276876) q[0];
rz(0.86979903) q[2];
sx q[2];
rz(-2.3491078) q[2];
sx q[2];
rz(0.078951702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15236552) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(2.3163124) q[1];
rz(0.18685347) q[3];
sx q[3];
rz(-2.0019751) q[3];
sx q[3];
rz(-0.29160515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(2.4287756) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7889439) q[0];
sx q[0];
rz(-1.7907356) q[0];
sx q[0];
rz(-0.21531944) q[0];
x q[1];
rz(-1.1586645) q[2];
sx q[2];
rz(-2.7325028) q[2];
sx q[2];
rz(1.8754174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7506999) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(0.23551029) q[1];
rz(-pi) q[2];
rz(2.0531822) q[3];
sx q[3];
rz(-1.6420806) q[3];
sx q[3];
rz(2.7977668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5793844) q[2];
sx q[2];
rz(0.66398579) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5749213) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(0.57408875) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-0.92837292) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1855477) q[0];
sx q[0];
rz(-1.569824) q[0];
sx q[0];
rz(2.2768343) q[0];
x q[1];
rz(1.9119772) q[2];
sx q[2];
rz(-0.8114292) q[2];
sx q[2];
rz(2.2974643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0005562) q[1];
sx q[1];
rz(-2.6489241) q[1];
sx q[1];
rz(-2.1129184) q[1];
rz(-pi) q[2];
rz(-2.6392322) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(-0.82752284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.1553923) q[2];
rz(-0.14285764) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(2.3891881) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518406) q[0];
sx q[0];
rz(-1.4983488) q[0];
sx q[0];
rz(-2.3206582) q[0];
x q[1];
rz(-0.98624595) q[2];
sx q[2];
rz(-2.2918321) q[2];
sx q[2];
rz(-2.4436488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4068027) q[1];
sx q[1];
rz(-0.64195913) q[1];
sx q[1];
rz(1.2924679) q[1];
rz(-pi) q[2];
rz(2.1390901) q[3];
sx q[3];
rz(-1.2270524) q[3];
sx q[3];
rz(-2.7078201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(-1.0460098) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9933269) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(2.6519725) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(1.6061868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.673296) q[0];
sx q[0];
rz(-1.522994) q[0];
sx q[0];
rz(0.099138069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3697482) q[2];
sx q[2];
rz(-2.2945885) q[2];
sx q[2];
rz(0.48230241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3121351) q[1];
sx q[1];
rz(-1.8541938) q[1];
sx q[1];
rz(-2.9764922) q[1];
x q[2];
rz(0.62932265) q[3];
sx q[3];
rz(-1.704708) q[3];
sx q[3];
rz(2.8117992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(1.2699132) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.3659182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89643909) q[0];
sx q[0];
rz(-1.23587) q[0];
sx q[0];
rz(2.4048637) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8132509) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(-2.6048425) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1508031) q[1];
sx q[1];
rz(-1.3152796) q[1];
sx q[1];
rz(1.7533416) q[1];
rz(-1.4468206) q[3];
sx q[3];
rz(-2.0967391) q[3];
sx q[3];
rz(2.7391609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0718677) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(2.243637) q[2];
rz(2.9979624) q[3];
sx q[3];
rz(-1.8840021) q[3];
sx q[3];
rz(-2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-0.32178497) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(2.5245573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8480393) q[0];
sx q[0];
rz(-0.32395054) q[0];
sx q[0];
rz(2.241961) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7865137) q[2];
sx q[2];
rz(-1.5383913) q[2];
sx q[2];
rz(2.0916794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1083793) q[1];
sx q[1];
rz(-1.2011659) q[1];
sx q[1];
rz(-2.0481678) q[1];
rz(0.7493895) q[3];
sx q[3];
rz(-1.6036311) q[3];
sx q[3];
rz(-2.4421305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(1.6581416) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9634919) q[0];
sx q[0];
rz(-3.0577457) q[0];
sx q[0];
rz(-1.447178) q[0];
x q[1];
rz(-2.9786417) q[2];
sx q[2];
rz(-2.5451247) q[2];
sx q[2];
rz(-0.60566723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9121453) q[1];
sx q[1];
rz(-1.7143237) q[1];
sx q[1];
rz(-1.9043282) q[1];
x q[2];
rz(0.55012196) q[3];
sx q[3];
rz(-1.9024444) q[3];
sx q[3];
rz(1.630993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1411529) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(1.2443939) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.72220951) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(1.3409021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6751624) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(1.947233) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4341899) q[2];
sx q[2];
rz(-2.280526) q[2];
sx q[2];
rz(-3.1106069) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(0.96149573) q[1];
rz(2.9756536) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(-0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.6293388) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.6188999) q[2];
rz(-0.86087888) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(-1.6265097) q[2];
sx q[2];
rz(-0.40691661) q[2];
sx q[2];
rz(-0.42452068) q[2];
rz(2.7545746) q[3];
sx q[3];
rz(-1.1292463) q[3];
sx q[3];
rz(-2.675727) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];