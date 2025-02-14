OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.926446) q[0];
sx q[0];
rz(-0.67987052) q[0];
sx q[0];
rz(3.0409066) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(-2.1159621) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2446063) q[0];
sx q[0];
rz(-1.375544) q[0];
sx q[0];
rz(-0.050873916) q[0];
x q[1];
rz(2.8335613) q[2];
sx q[2];
rz(-1.5218166) q[2];
sx q[2];
rz(0.9020976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.034596171) q[1];
sx q[1];
rz(-1.6507848) q[1];
sx q[1];
rz(-0.75350113) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64798665) q[3];
sx q[3];
rz(-0.77666908) q[3];
sx q[3];
rz(2.4512993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2166298) q[2];
sx q[2];
rz(-1.3684042) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(1.9803068) q[3];
sx q[3];
rz(-2.1372644) q[3];
sx q[3];
rz(-1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0733114) q[0];
sx q[0];
rz(-2.735205) q[0];
sx q[0];
rz(2.708129) q[0];
rz(2.0794226) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(-2.4210222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90209157) q[0];
sx q[0];
rz(-1.4441617) q[0];
sx q[0];
rz(0.037221639) q[0];
rz(1.9754875) q[2];
sx q[2];
rz(-2.5705159) q[2];
sx q[2];
rz(0.80246682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1432508) q[1];
sx q[1];
rz(-1.2187119) q[1];
sx q[1];
rz(-1.7290753) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0508762) q[3];
sx q[3];
rz(-0.54819104) q[3];
sx q[3];
rz(-0.91778008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34708193) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(-2.1573055) q[2];
rz(2.2847564) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(-1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03265753) q[0];
sx q[0];
rz(-2.8553243) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(3.1318829) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3056639) q[0];
sx q[0];
rz(-1.9171414) q[0];
sx q[0];
rz(0.54327528) q[0];
x q[1];
rz(2.9417324) q[2];
sx q[2];
rz(-1.4132573) q[2];
sx q[2];
rz(0.021856088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.328542) q[1];
sx q[1];
rz(-2.790741) q[1];
sx q[1];
rz(2.0787129) q[1];
rz(-0.023223485) q[3];
sx q[3];
rz(-1.6722213) q[3];
sx q[3];
rz(1.9254799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0297086) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(0.65579826) q[2];
rz(-0.2119952) q[3];
sx q[3];
rz(-0.63181221) q[3];
sx q[3];
rz(-1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(-1.5546881) q[0];
rz(0.2298062) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.302964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6617216) q[0];
sx q[0];
rz(-2.4085143) q[0];
sx q[0];
rz(0.69723155) q[0];
rz(-0.83452819) q[2];
sx q[2];
rz(-0.38943651) q[2];
sx q[2];
rz(1.8846041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6524383) q[1];
sx q[1];
rz(-1.592318) q[1];
sx q[1];
rz(-2.6179594) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0423562) q[3];
sx q[3];
rz(-2.5109249) q[3];
sx q[3];
rz(1.809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.135123) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(0.50372493) q[2];
rz(-1.5375562) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(1.5639719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85394323) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(2.8751539) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(2.0727167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7212413) q[0];
sx q[0];
rz(-1.9801894) q[0];
sx q[0];
rz(2.4024525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8891625) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(-2.3311262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0068897) q[1];
sx q[1];
rz(-1.1313032) q[1];
sx q[1];
rz(-1.2864) q[1];
rz(-pi) q[2];
rz(1.8183299) q[3];
sx q[3];
rz(-1.8283852) q[3];
sx q[3];
rz(3.0158693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2499007) q[2];
sx q[2];
rz(-2.4155858) q[2];
sx q[2];
rz(0.1740087) q[2];
rz(1.2627259) q[3];
sx q[3];
rz(-1.6014674) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0083017666) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(-2.2513576) q[0];
rz(1.4324191) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(-0.11996809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9667119) q[0];
sx q[0];
rz(-2.565756) q[0];
sx q[0];
rz(-0.26782803) q[0];
rz(-1.1402848) q[2];
sx q[2];
rz(-1.8906279) q[2];
sx q[2];
rz(-2.8240273) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6585582) q[1];
sx q[1];
rz(-2.2081828) q[1];
sx q[1];
rz(-1.263522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8268485) q[3];
sx q[3];
rz(-2.5786244) q[3];
sx q[3];
rz(1.1551577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.031781901) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(2.7500582) q[2];
rz(1.3854965) q[3];
sx q[3];
rz(-0.23844312) q[3];
sx q[3];
rz(-0.63867205) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077797) q[0];
sx q[0];
rz(-2.5901828) q[0];
sx q[0];
rz(-1.8361924) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69687122) q[0];
sx q[0];
rz(-1.5550647) q[0];
sx q[0];
rz(-3.12495) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6487892) q[2];
sx q[2];
rz(-1.7597919) q[2];
sx q[2];
rz(-2.1616621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0266155) q[1];
sx q[1];
rz(-1.3040286) q[1];
sx q[1];
rz(1.2342427) q[1];
rz(-pi) q[2];
rz(-1.9402294) q[3];
sx q[3];
rz(-1.2681586) q[3];
sx q[3];
rz(-1.7517881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(-3.0090561) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(-2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9082311) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(0.37286266) q[0];
rz(0.33509675) q[1];
sx q[1];
rz(-1.2058039) q[1];
sx q[1];
rz(-2.0705409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.459254) q[0];
sx q[0];
rz(-1.5821348) q[0];
sx q[0];
rz(-3.0724694) q[0];
x q[1];
rz(2.5899289) q[2];
sx q[2];
rz(-1.7529738) q[2];
sx q[2];
rz(-0.38924402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0275299) q[1];
sx q[1];
rz(-0.34110554) q[1];
sx q[1];
rz(0.76082629) q[1];
x q[2];
rz(3.0464977) q[3];
sx q[3];
rz(-1.250306) q[3];
sx q[3];
rz(1.8959498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7143453) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(1.0672807) q[2];
rz(2.769477) q[3];
sx q[3];
rz(-2.1477063) q[3];
sx q[3];
rz(0.20017643) q[3];
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
rz(-pi) q[3];
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
rz(-1.9099971) q[0];
sx q[0];
rz(-2.5439926) q[0];
sx q[0];
rz(-1.4071314) q[0];
rz(-2.097997) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(2.6503906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1192899) q[0];
sx q[0];
rz(-0.91141846) q[0];
sx q[0];
rz(-2.7727866) q[0];
rz(2.3517026) q[2];
sx q[2];
rz(-2.9040376) q[2];
sx q[2];
rz(1.6611163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8206827) q[1];
sx q[1];
rz(-1.2315799) q[1];
sx q[1];
rz(-1.7103819) q[1];
rz(-pi) q[2];
rz(0.43213958) q[3];
sx q[3];
rz(-0.41101217) q[3];
sx q[3];
rz(-0.43288818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6970984) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(2.3704884) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412909) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(1.6218761) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.8404574) q[1];
sx q[1];
rz(-0.69449743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3169169) q[0];
sx q[0];
rz(-1.7986713) q[0];
sx q[0];
rz(0.014046305) q[0];
x q[1];
rz(2.6171309) q[2];
sx q[2];
rz(-1.6860644) q[2];
sx q[2];
rz(-1.0712717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2934781) q[1];
sx q[1];
rz(-1.2367931) q[1];
sx q[1];
rz(3.1298742) q[1];
x q[2];
rz(1.1841838) q[3];
sx q[3];
rz(-0.85672934) q[3];
sx q[3];
rz(2.3165645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0684315) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(2.8237776) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(-2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(1.9152676) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(2.4438159) q[2];
sx q[2];
rz(-1.991376) q[2];
sx q[2];
rz(-2.1817653) q[2];
rz(2.0967284) q[3];
sx q[3];
rz(-1.5689701) q[3];
sx q[3];
rz(0.29472385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
