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
rz(1.7841568) q[0];
sx q[0];
rz(-0.62499243) q[0];
sx q[0];
rz(-1.619119) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(3.7781236) q[1];
sx q[1];
rz(14.349024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5769437) q[0];
sx q[0];
rz(-1.859382) q[0];
sx q[0];
rz(2.1493069) q[0];
x q[1];
rz(2.236249) q[2];
sx q[2];
rz(-0.45768379) q[2];
sx q[2];
rz(-2.8913218) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57558838) q[1];
sx q[1];
rz(-1.8270565) q[1];
sx q[1];
rz(-2.9017067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3728849) q[3];
sx q[3];
rz(-2.6555914) q[3];
sx q[3];
rz(2.3305691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4437359) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(3.0016628) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(-0.073277624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(2.7781558) q[0];
rz(-2.4711171) q[1];
sx q[1];
rz(-1.3251708) q[1];
sx q[1];
rz(-1.0796116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19980656) q[0];
sx q[0];
rz(-0.58871709) q[0];
sx q[0];
rz(2.0090282) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.011376) q[2];
sx q[2];
rz(-1.0417582) q[2];
sx q[2];
rz(2.0646273) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9448589) q[1];
sx q[1];
rz(-0.24538876) q[1];
sx q[1];
rz(1.3887482) q[1];
rz(-0.58066253) q[3];
sx q[3];
rz(-1.7911909) q[3];
sx q[3];
rz(3.0739258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4662027) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(-2.7030763) q[2];
rz(1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587875) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(-1.8939182) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(2.8819328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369401) q[0];
sx q[0];
rz(-0.90422179) q[0];
sx q[0];
rz(-2.7949692) q[0];
rz(-pi) q[1];
rz(-2.3724062) q[2];
sx q[2];
rz(-2.63317) q[2];
sx q[2];
rz(1.2505797) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34777113) q[1];
sx q[1];
rz(-2.035537) q[1];
sx q[1];
rz(-0.79885599) q[1];
rz(-pi) q[2];
rz(0.91783701) q[3];
sx q[3];
rz(-1.3363095) q[3];
sx q[3];
rz(-3.0378436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20716771) q[2];
sx q[2];
rz(-1.5639037) q[2];
sx q[2];
rz(2.3210607) q[2];
rz(-0.82646838) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(1.8291399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25927037) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(-0.93233863) q[0];
rz(1.3177634) q[1];
sx q[1];
rz(-1.1612786) q[1];
sx q[1];
rz(-0.12277776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6500191) q[0];
sx q[0];
rz(-2.5959925) q[0];
sx q[0];
rz(0.16262098) q[0];
x q[1];
rz(-0.81715314) q[2];
sx q[2];
rz(-1.3434402) q[2];
sx q[2];
rz(1.0447234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40298324) q[1];
sx q[1];
rz(-1.5697641) q[1];
sx q[1];
rz(2.3139075) q[1];
x q[2];
rz(1.8639542) q[3];
sx q[3];
rz(-1.0256301) q[3];
sx q[3];
rz(2.4923862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54470283) q[2];
sx q[2];
rz(-0.92146102) q[2];
sx q[2];
rz(2.5370562) q[2];
rz(-2.4321411) q[3];
sx q[3];
rz(-2.3716898) q[3];
sx q[3];
rz(2.3179222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17968793) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(-1.3667579) q[0];
rz(-2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(2.6008115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7779988) q[0];
sx q[0];
rz(-2.8907597) q[0];
sx q[0];
rz(0.59441113) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3983058) q[2];
sx q[2];
rz(-1.0096482) q[2];
sx q[2];
rz(0.66711344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4879228) q[1];
sx q[1];
rz(-1.2752295) q[1];
sx q[1];
rz(-1.3699156) q[1];
x q[2];
rz(2.861372) q[3];
sx q[3];
rz(-2.5570405) q[3];
sx q[3];
rz(1.9909797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6775386) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(-0.82826725) q[2];
rz(-1.4520491) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(0.90788466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.117711) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(1.1718132) q[0];
rz(2.4614629) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(-0.5562869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2883462) q[0];
sx q[0];
rz(-2.6814406) q[0];
sx q[0];
rz(-0.96921885) q[0];
rz(-pi) q[1];
rz(1.589619) q[2];
sx q[2];
rz(-0.14655098) q[2];
sx q[2];
rz(2.3702247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7032675) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(3.0238815) q[1];
rz(-pi) q[2];
rz(-2.049796) q[3];
sx q[3];
rz(-0.99078876) q[3];
sx q[3];
rz(-0.73778462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3351626) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(-0.20509091) q[2];
rz(2.3388376) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(0.55454379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270808) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(-2.9141973) q[0];
rz(-2.3511476) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(2.3048185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6444708) q[0];
sx q[0];
rz(-2.5724412) q[0];
sx q[0];
rz(2.7272124) q[0];
rz(-pi) q[1];
rz(1.2143986) q[2];
sx q[2];
rz(-1.0123535) q[2];
sx q[2];
rz(1.7624378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6001125) q[1];
sx q[1];
rz(-1.8842234) q[1];
sx q[1];
rz(0.42509834) q[1];
rz(0.76117875) q[3];
sx q[3];
rz(-1.8463148) q[3];
sx q[3];
rz(2.8305284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5493912) q[2];
sx q[2];
rz(-1.4791919) q[2];
sx q[2];
rz(0.57360348) q[2];
rz(0.32522374) q[3];
sx q[3];
rz(-0.89541382) q[3];
sx q[3];
rz(-0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-0.17806299) q[0];
sx q[0];
rz(0.67894116) q[0];
rz(-3.0067054) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(1.4237684) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2401827) q[0];
sx q[0];
rz(-1.115832) q[0];
sx q[0];
rz(-0.97413537) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5871954) q[2];
sx q[2];
rz(-0.94292484) q[2];
sx q[2];
rz(0.45931057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9805124) q[1];
sx q[1];
rz(-2.3337165) q[1];
sx q[1];
rz(0.041461583) q[1];
rz(-pi) q[2];
rz(-1.1066439) q[3];
sx q[3];
rz(-1.3106133) q[3];
sx q[3];
rz(-1.2342208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5941102) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(-2.9098848) q[2];
rz(1.0101275) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(0.89099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6758839) q[0];
sx q[0];
rz(-2.3864585) q[0];
sx q[0];
rz(0.51374197) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(1.5076216) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9378273) q[0];
sx q[0];
rz(-1.7047593) q[0];
sx q[0];
rz(-1.521827) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7319674) q[2];
sx q[2];
rz(-1.8212166) q[2];
sx q[2];
rz(2.7572244) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3097174) q[1];
sx q[1];
rz(-0.17942218) q[1];
sx q[1];
rz(2.37876) q[1];
x q[2];
rz(-1.6674897) q[3];
sx q[3];
rz(-2.8887199) q[3];
sx q[3];
rz(-2.0677572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12019176) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(-2.8821442) q[2];
rz(-1.1546968) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(1.8416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65086377) q[0];
sx q[0];
rz(-1.4800973) q[0];
sx q[0];
rz(-1.65253) q[0];
rz(2.5015855) q[1];
sx q[1];
rz(-2.044544) q[1];
sx q[1];
rz(-0.99001137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418808) q[0];
sx q[0];
rz(-0.16946259) q[0];
sx q[0];
rz(-2.5925596) q[0];
x q[1];
rz(-2.3732164) q[2];
sx q[2];
rz(-1.2207165) q[2];
sx q[2];
rz(-2.1742252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23259904) q[1];
sx q[1];
rz(-0.4505583) q[1];
sx q[1];
rz(1.1934936) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1166273) q[3];
sx q[3];
rz(-1.7933266) q[3];
sx q[3];
rz(-0.74786598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9330357) q[2];
sx q[2];
rz(-1.2988337) q[2];
sx q[2];
rz(0.2253069) q[2];
rz(2.2202668) q[3];
sx q[3];
rz(-2.7511629) q[3];
sx q[3];
rz(-3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13023547) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(-1.8232952) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(3.0443424) q[2];
sx q[2];
rz(-1.3122059) q[2];
sx q[2];
rz(-2.9065804) q[2];
rz(-2.1193567) q[3];
sx q[3];
rz(-1.7283741) q[3];
sx q[3];
rz(-1.7909539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
