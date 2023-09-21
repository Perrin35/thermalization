OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0041381) q[0];
sx q[0];
rz(-1.9507427) q[0];
sx q[0];
rz(-0.11560346) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9723188) q[2];
sx q[2];
rz(-1.2520773) q[2];
sx q[2];
rz(-0.67414325) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5296386) q[1];
sx q[1];
rz(-2.0679592) q[1];
sx q[1];
rz(-2.0938718) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8205809) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(0.111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(-1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(2.3056727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3859235) q[0];
sx q[0];
rz(-0.14005157) q[0];
sx q[0];
rz(1.7007909) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0216653) q[2];
sx q[2];
rz(-1.2163391) q[2];
sx q[2];
rz(-0.82565386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4513431) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(0.33957014) q[1];
x q[2];
rz(-0.73987506) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(1.5996931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(2.704481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8374098) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(-3.0512179) q[0];
rz(-1.978546) q[2];
sx q[2];
rz(-1.7951269) q[2];
sx q[2];
rz(0.11220223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6091842) q[1];
sx q[1];
rz(-2.6951365) q[1];
sx q[1];
rz(-2.8162454) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4594853) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(1.5745844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(-1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(0.19128004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080973074) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(3.0540375) q[0];
rz(-pi) q[1];
rz(-1.8727826) q[2];
sx q[2];
rz(-2.3846845) q[2];
sx q[2];
rz(2.7243171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0335238) q[1];
sx q[1];
rz(-2.3571157) q[1];
sx q[1];
rz(-0.14524059) q[1];
rz(-pi) q[2];
rz(2.6077765) q[3];
sx q[3];
rz(-1.8752408) q[3];
sx q[3];
rz(-2.387405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33870944) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-2.5849735) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(-0.97250485) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(1.6385965) q[0];
x q[1];
rz(2.3210578) q[2];
sx q[2];
rz(-1.8168601) q[2];
sx q[2];
rz(-0.29124242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16469615) q[1];
sx q[1];
rz(-1.7105192) q[1];
sx q[1];
rz(-1.3752027) q[1];
x q[2];
rz(-2.186741) q[3];
sx q[3];
rz(-1.6007489) q[3];
sx q[3];
rz(2.0314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(2.1877066) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4218629) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(0.37429601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930646) q[0];
sx q[0];
rz(-2.6120798) q[0];
sx q[0];
rz(-0.66324309) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6475545) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(-1.5649232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(-2.5724263) q[1];
x q[2];
rz(0.88126392) q[3];
sx q[3];
rz(-1.0242108) q[3];
sx q[3];
rz(-1.5342086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(1.8315171) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3416672) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-0.13473405) q[0];
x q[1];
rz(-2.4285165) q[2];
sx q[2];
rz(-1.3165858) q[2];
sx q[2];
rz(1.6377246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3771364) q[1];
sx q[1];
rz(-2.7546282) q[1];
sx q[1];
rz(-0.52486921) q[1];
rz(-pi) q[2];
rz(1.5433611) q[3];
sx q[3];
rz(-1.0155639) q[3];
sx q[3];
rz(-2.1813986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2157796) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(2.8619134) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(2.2917152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4445164) q[0];
sx q[0];
rz(-2.8651617) q[0];
sx q[0];
rz(-2.1806549) q[0];
rz(-pi) q[1];
rz(3.0476961) q[2];
sx q[2];
rz(-0.88985032) q[2];
sx q[2];
rz(-2.4927793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0975768) q[1];
sx q[1];
rz(-2.3313287) q[1];
sx q[1];
rz(0.047659831) q[1];
rz(-pi) q[2];
rz(-2.0073118) q[3];
sx q[3];
rz(-0.28537649) q[3];
sx q[3];
rz(-2.7989822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.3809416) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(-1.2776432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.02361) q[0];
sx q[0];
rz(-1.3408459) q[0];
sx q[0];
rz(0.50595565) q[0];
rz(-pi) q[1];
rz(-0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(0.28011766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9459878) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(0.040954879) q[1];
rz(-pi) q[2];
rz(-0.38655917) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(-2.3545635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-0.35783106) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-1.0740124) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-0.11225587) q[0];
rz(0.90011251) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(-2.9311438) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68073273) q[0];
sx q[0];
rz(-1.0675149) q[0];
sx q[0];
rz(-1.3637278) q[0];
rz(-pi) q[1];
rz(-1.9333282) q[2];
sx q[2];
rz(-2.7047485) q[2];
sx q[2];
rz(-0.7561965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12293832) q[1];
sx q[1];
rz(-1.2631053) q[1];
sx q[1];
rz(1.3535181) q[1];
rz(-pi) q[2];
rz(1.8080766) q[3];
sx q[3];
rz(-1.3271601) q[3];
sx q[3];
rz(-2.7436649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18008733) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(-2.9755637) q[2];
sx q[2];
rz(-2.0353073) q[2];
sx q[2];
rz(-1.4370949) q[2];
rz(-1.9498701) q[3];
sx q[3];
rz(-0.76615292) q[3];
sx q[3];
rz(-2.5828008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];