OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(0.83067218) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2031189) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(-1.9012326) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4029998) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(1.5473168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9484529) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(2.302041) q[1];
x q[2];
rz(1.9740231) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(-2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3347496) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1268483) q[2];
sx q[2];
rz(-1.7146646) q[2];
sx q[2];
rz(0.680188) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3152299) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(-2.1748494) q[1];
rz(-pi) q[2];
rz(-1.4665524) q[3];
sx q[3];
rz(-0.6904656) q[3];
sx q[3];
rz(-3.1363827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49845093) q[0];
sx q[0];
rz(-2.8006449) q[0];
sx q[0];
rz(-1.5091512) q[0];
rz(-pi) q[1];
rz(-1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(2.0470326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87103802) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(0.87267455) q[1];
rz(-pi) q[2];
rz(0.4226513) q[3];
sx q[3];
rz(-1.7461516) q[3];
sx q[3];
rz(-2.5734176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-2.9761956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9497313) q[0];
sx q[0];
rz(-0.22338578) q[0];
sx q[0];
rz(2.1641157) q[0];
x q[1];
rz(2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-0.93726678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0927825) q[1];
sx q[1];
rz(-2.4895992) q[1];
sx q[1];
rz(-0.55744967) q[1];
rz(-pi) q[2];
rz(-1.0159675) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-2.9096471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(2.288726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69552341) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-2.7930789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2105337) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(1.7935497) q[1];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4076685) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(2.2785447) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73156725) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(-2.3716795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62085405) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(1.1362856) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44316767) q[3];
sx q[3];
rz(-0.89649761) q[3];
sx q[3];
rz(-2.5600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543024) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.2994231) q[0];
rz(-pi) q[1];
rz(-1.5067528) q[2];
sx q[2];
rz(-1.8539068) q[2];
sx q[2];
rz(1.7664906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1544513) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229066) q[0];
sx q[0];
rz(-1.6081928) q[0];
sx q[0];
rz(-0.54625578) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9756873) q[2];
sx q[2];
rz(-0.84120175) q[2];
sx q[2];
rz(-1.0798432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(-0.016896292) q[1];
x q[2];
rz(0.47655388) q[3];
sx q[3];
rz(-1.8503975) q[3];
sx q[3];
rz(2.9001146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-2.6064176) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(-0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(-2.2156782) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(2.9979343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70897734) q[2];
sx q[2];
rz(-0.93009863) q[2];
sx q[2];
rz(0.35352732) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.029286413) q[1];
sx q[1];
rz(-0.46488133) q[1];
sx q[1];
rz(-1.2390562) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52243201) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-1.013247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(2.9985715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(1.313414) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71868616) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(-2.7194517) q[1];
x q[2];
rz(0.69443955) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(-1.6758067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(1.1258833) q[2];
sx q[2];
rz(-1.6402259) q[2];
sx q[2];
rz(0.31625407) q[2];
rz(3.0782386) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
