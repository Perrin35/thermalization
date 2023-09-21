OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8496348) q[0];
sx q[0];
rz(-0.44009527) q[0];
sx q[0];
rz(0.13719288) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(-1.7383716) q[1];
sx q[1];
rz(0.52991968) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3903025) q[0];
sx q[0];
rz(-1.4634702) q[0];
sx q[0];
rz(-1.9530549) q[0];
x q[1];
rz(2.2719703) q[2];
sx q[2];
rz(-0.50719559) q[2];
sx q[2];
rz(-1.6091572) q[2];
rz(-pi) q[3];
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
rz(1.3210117) q[3];
sx q[3];
rz(-2.3134109) q[3];
sx q[3];
rz(0.111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(1.9477828) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(-2.3056727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7556691) q[0];
sx q[0];
rz(-3.0015411) q[0];
sx q[0];
rz(-1.4408017) q[0];
rz(-pi) q[1];
rz(2.732227) q[2];
sx q[2];
rz(-1.0592807) q[2];
sx q[2];
rz(-2.6057233) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.086120124) q[1];
sx q[1];
rz(-1.8430084) q[1];
sx q[1];
rz(-0.91113669) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73987506) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-1.0536449) q[0];
rz(1.9127649) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(-2.704481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99840435) q[0];
sx q[0];
rz(-0.30448738) q[0];
sx q[0];
rz(-1.2782774) q[0];
x q[1];
rz(1.978546) q[2];
sx q[2];
rz(-1.3464658) q[2];
sx q[2];
rz(0.11220223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8902953) q[1];
sx q[1];
rz(-1.9922868) q[1];
sx q[1];
rz(-1.7226268) q[1];
x q[2];
rz(-1.0760355) q[3];
sx q[3];
rz(-1.6238188) q[3];
sx q[3];
rz(0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-2.9503126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6879038) q[0];
sx q[0];
rz(-2.9691302) q[0];
sx q[0];
rz(2.0989291) q[0];
x q[1];
rz(-1.26881) q[2];
sx q[2];
rz(-0.75690818) q[2];
sx q[2];
rz(2.7243171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8297255) q[1];
sx q[1];
rz(-0.79677454) q[1];
sx q[1];
rz(1.7142678) q[1];
x q[2];
rz(2.6077765) q[3];
sx q[3];
rz(-1.2663519) q[3];
sx q[3];
rz(2.387405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(2.130924) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33870944) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(2.5849735) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(1.6385965) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9237256) q[2];
sx q[2];
rz(-0.7820411) q[2];
sx q[2];
rz(1.023934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9768965) q[1];
sx q[1];
rz(-1.7105192) q[1];
sx q[1];
rz(-1.7663899) q[1];
x q[2];
rz(3.1049018) q[3];
sx q[3];
rz(-2.1864236) q[3];
sx q[3];
rz(2.7021367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(-1.3657773) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.4259889) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(0.37429601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22898856) q[0];
sx q[0];
rz(-1.8870263) q[0];
sx q[0];
rz(2.7094748) q[0];
rz(-pi) q[1];
rz(-0.46005581) q[2];
sx q[2];
rz(-0.54121491) q[2];
sx q[2];
rz(2.7341026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7172076) q[1];
sx q[1];
rz(-1.3076412) q[1];
sx q[1];
rz(2.5724263) q[1];
rz(-2.2603287) q[3];
sx q[3];
rz(-2.1173819) q[3];
sx q[3];
rz(1.5342086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(-0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-2.5206101) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.8315171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.037127) q[0];
sx q[0];
rz(-1.5195527) q[0];
sx q[0];
rz(2.7535704) q[0];
x q[1];
rz(1.2398948) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(-2.9943525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.29855) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(-0.33903867) q[1];
rz(-1.5982315) q[3];
sx q[3];
rz(-1.0155639) q[3];
sx q[3];
rz(-2.1813986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.6581992) q[2];
rz(2.8619134) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-0.057597615) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(-2.9220007) q[0];
rz(-2.638468) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(-2.2917152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69707623) q[0];
sx q[0];
rz(-2.8651617) q[0];
sx q[0];
rz(-2.1806549) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6860028) q[2];
sx q[2];
rz(-0.68636471) q[2];
sx q[2];
rz(-2.3442868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1666607) q[1];
sx q[1];
rz(-2.3798675) q[1];
sx q[1];
rz(1.5207661) q[1];
rz(3.0181846) q[3];
sx q[3];
rz(-1.3128237) q[3];
sx q[3];
rz(2.3464399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(-2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5091771) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(2.6889154) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(1.2776432) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15690878) q[0];
sx q[0];
rz(-2.589993) q[0];
sx q[0];
rz(0.44996913) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56565506) q[2];
sx q[2];
rz(-2.3662162) q[2];
sx q[2];
rz(2.2765809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
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
x q[2];
rz(-0.38655917) q[3];
sx q[3];
rz(-0.87214008) q[3];
sx q[3];
rz(2.3545635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8032802) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.91838592) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-0.11225587) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(0.21044883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716547) q[0];
sx q[0];
rz(-0.5408113) q[0];
sx q[0];
rz(-0.35738118) q[0];
x q[1];
rz(2.9774882) q[2];
sx q[2];
rz(-1.1640942) q[2];
sx q[2];
rz(1.9889529) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0186543) q[1];
sx q[1];
rz(-1.8784874) q[1];
sx q[1];
rz(-1.7880746) q[1];
rz(-pi) q[2];
rz(0.75746234) q[3];
sx q[3];
rz(-0.33843455) q[3];
sx q[3];
rz(0.38871845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(-1.2522092) q[2];
sx q[2];
rz(-0.49124419) q[2];
sx q[2];
rz(-1.0791525) q[2];
rz(2.7995085) q[3];
sx q[3];
rz(-2.2707006) q[3];
sx q[3];
rz(-2.0778098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];