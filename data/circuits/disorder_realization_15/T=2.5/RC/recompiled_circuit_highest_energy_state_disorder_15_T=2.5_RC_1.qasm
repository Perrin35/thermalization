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
rz(-2.5315142) q[0];
sx q[0];
rz(-0.22992034) q[0];
sx q[0];
rz(0.71969405) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(0.12535867) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7490112) q[0];
sx q[0];
rz(-1.9960989) q[0];
sx q[0];
rz(-1.1293189) q[0];
rz(1.4761476) q[2];
sx q[2];
rz(-1.2834335) q[2];
sx q[2];
rz(-0.3269302) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64588378) q[1];
sx q[1];
rz(-1.570374) q[1];
sx q[1];
rz(-1.5721815) q[1];
x q[2];
rz(-0.35644021) q[3];
sx q[3];
rz(-2.2619777) q[3];
sx q[3];
rz(2.3158642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7980935) q[2];
sx q[2];
rz(-2.7332879) q[2];
sx q[2];
rz(-0.84585345) q[2];
rz(-2.3503303) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(0.066789269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(-0.043664232) q[0];
rz(1.5664258) q[1];
sx q[1];
rz(-1.7685726) q[1];
sx q[1];
rz(1.498819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18066809) q[0];
sx q[0];
rz(-2.5378413) q[0];
sx q[0];
rz(-1.5779499) q[0];
x q[1];
rz(-3.1286521) q[2];
sx q[2];
rz(-0.99943752) q[2];
sx q[2];
rz(0.017627942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0017569204) q[1];
sx q[1];
rz(-0.082232177) q[1];
sx q[1];
rz(-2.7539192) q[1];
rz(-pi) q[2];
x q[2];
rz(0.033229251) q[3];
sx q[3];
rz(-1.9150315) q[3];
sx q[3];
rz(-2.0118464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.046039) q[2];
sx q[2];
rz(-0.150103) q[2];
sx q[2];
rz(-2.6138439) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(-1.2214448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7432778) q[0];
sx q[0];
rz(-0.97284955) q[0];
sx q[0];
rz(-1.1463746) q[0];
rz(1.3829117) q[1];
sx q[1];
rz(-0.29255602) q[1];
sx q[1];
rz(3.0376099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706174) q[0];
sx q[0];
rz(-1.3461324) q[0];
sx q[0];
rz(-1.490834) q[0];
rz(-0.017113233) q[2];
sx q[2];
rz(-1.2984973) q[2];
sx q[2];
rz(2.2382617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3465991) q[1];
sx q[1];
rz(-1.775034) q[1];
sx q[1];
rz(-2.380382) q[1];
x q[2];
rz(-2.7894734) q[3];
sx q[3];
rz(-2.8330292) q[3];
sx q[3];
rz(-1.7087913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18342239) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(2.5799694) q[2];
rz(3.0567452) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(0.00076278846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7247923) q[0];
sx q[0];
rz(-2.8775207) q[0];
sx q[0];
rz(-2.6208139) q[0];
rz(2.9837823) q[1];
sx q[1];
rz(-2.4746555) q[1];
sx q[1];
rz(0.075798362) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2139839) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(1.7068205) q[0];
x q[1];
rz(0.00064223991) q[2];
sx q[2];
rz(-1.5721605) q[2];
sx q[2];
rz(1.431992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6356018) q[1];
sx q[1];
rz(-1.4293993) q[1];
sx q[1];
rz(1.1902203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10785477) q[3];
sx q[3];
rz(-1.1198992) q[3];
sx q[3];
rz(-1.0552561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51125222) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(1.7208257) q[2];
rz(-0.00682791) q[3];
sx q[3];
rz(-0.029191645) q[3];
sx q[3];
rz(-1.487287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0066978) q[0];
sx q[0];
rz(-2.9338574) q[0];
sx q[0];
rz(-0.2625221) q[0];
rz(0.94995704) q[1];
sx q[1];
rz(-3.0634395) q[1];
sx q[1];
rz(0.21771678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6852606) q[0];
sx q[0];
rz(-1.1781516) q[0];
sx q[0];
rz(-0.13360191) q[0];
rz(-pi) q[1];
rz(-2.3764059) q[2];
sx q[2];
rz(-3.0209241) q[2];
sx q[2];
rz(2.3007002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1982983) q[1];
sx q[1];
rz(-1.5729331) q[1];
sx q[1];
rz(1.398477) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0721139) q[3];
sx q[3];
rz(-1.0437878) q[3];
sx q[3];
rz(-1.219092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21506423) q[2];
sx q[2];
rz(-0.0096409163) q[2];
sx q[2];
rz(0.24791524) q[2];
rz(-1.3748112) q[3];
sx q[3];
rz(-0.050516613) q[3];
sx q[3];
rz(1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1011825) q[0];
sx q[0];
rz(-1.269416) q[0];
sx q[0];
rz(-0.61224473) q[0];
rz(-0.170389) q[1];
sx q[1];
rz(-3.0590765) q[1];
sx q[1];
rz(1.4917397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1338324) q[0];
sx q[0];
rz(-0.88307021) q[0];
sx q[0];
rz(-0.69761116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5666008) q[2];
sx q[2];
rz(-1.5853264) q[2];
sx q[2];
rz(-0.57455237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41497624) q[1];
sx q[1];
rz(-1.5041989) q[1];
sx q[1];
rz(2.8482751) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0759367) q[3];
sx q[3];
rz(-2.7759984) q[3];
sx q[3];
rz(-2.5557809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49543574) q[2];
sx q[2];
rz(-3.1313681) q[2];
sx q[2];
rz(0.23908991) q[2];
rz(0.80860364) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(-1.3330207) q[3];
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
rz(-pi/2) q[0];
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
rz(2.7072153) q[0];
sx q[0];
rz(-3.0476397) q[0];
sx q[0];
rz(0.77675003) q[0];
rz(0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.6687261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4368255) q[0];
sx q[0];
rz(-0.87427199) q[0];
sx q[0];
rz(-0.27946194) q[0];
rz(-0.02350925) q[2];
sx q[2];
rz(-1.5644367) q[2];
sx q[2];
rz(1.7451177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1030813) q[1];
sx q[1];
rz(-2.3902479) q[1];
sx q[1];
rz(2.9935212) q[1];
rz(-0.075447791) q[3];
sx q[3];
rz(-1.8561072) q[3];
sx q[3];
rz(-2.4782654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1482859) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(-0.41603184) q[2];
rz(-2.3038583) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(-1.4778888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763022) q[0];
sx q[0];
rz(-3.1017113) q[0];
sx q[0];
rz(-0.95529977) q[0];
rz(-1.2166066) q[1];
sx q[1];
rz(-0.33733264) q[1];
sx q[1];
rz(-2.7359656) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8389616) q[0];
sx q[0];
rz(-1.8895188) q[0];
sx q[0];
rz(-0.74619729) q[0];
rz(-0.50875278) q[2];
sx q[2];
rz(-0.18157427) q[2];
sx q[2];
rz(-1.5588829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9863524) q[1];
sx q[1];
rz(-1.494368) q[1];
sx q[1];
rz(-2.9806311) q[1];
x q[2];
rz(0.38281103) q[3];
sx q[3];
rz(-0.22914003) q[3];
sx q[3];
rz(-2.7349796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8568628) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(-2.3178318) q[2];
rz(-0.13122261) q[3];
sx q[3];
rz(-0.28990144) q[3];
sx q[3];
rz(0.32740617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.71309483) q[0];
sx q[0];
rz(-0.14398028) q[0];
sx q[0];
rz(2.4289828) q[0];
rz(-2.5932942) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(-0.20864329) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8241103) q[0];
sx q[0];
rz(-1.1718996) q[0];
sx q[0];
rz(-2.0092416) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1227112) q[2];
sx q[2];
rz(-1.5361726) q[2];
sx q[2];
rz(-2.0378626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.029034792) q[1];
sx q[1];
rz(-2.1264972) q[1];
sx q[1];
rz(-0.0056000756) q[1];
rz(-pi) q[2];
rz(-2.3712158) q[3];
sx q[3];
rz(-0.80157864) q[3];
sx q[3];
rz(1.0018536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.448552) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(2.9244259) q[2];
rz(2.7740357) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(2.1448081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17598584) q[0];
sx q[0];
rz(-0.088986926) q[0];
sx q[0];
rz(-0.17372818) q[0];
rz(-1.732775) q[1];
sx q[1];
rz(-1.4585739) q[1];
sx q[1];
rz(1.4558314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98783606) q[0];
sx q[0];
rz(-2.5964886) q[0];
sx q[0];
rz(-2.821652) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0020724) q[2];
sx q[2];
rz(-1.7762842) q[2];
sx q[2];
rz(0.71256283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10911726) q[1];
sx q[1];
rz(-0.87504234) q[1];
sx q[1];
rz(-1.7233707) q[1];
x q[2];
rz(-0.10723857) q[3];
sx q[3];
rz(-1.3235561) q[3];
sx q[3];
rz(0.040217248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(-1.3803587) q[2];
rz(-0.68584758) q[3];
sx q[3];
rz(-0.0018456056) q[3];
sx q[3];
rz(2.4628911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.3995517) q[0];
sx q[0];
rz(-2.4492332) q[0];
sx q[0];
rz(-1.5204182) q[0];
rz(1.5581268) q[1];
sx q[1];
rz(-1.5065267) q[1];
sx q[1];
rz(-2.9361257) q[1];
rz(-1.6995399) q[2];
sx q[2];
rz(-3.0151571) q[2];
sx q[2];
rz(-2.7966316) q[2];
rz(-3.0229983) q[3];
sx q[3];
rz(-1.8684602) q[3];
sx q[3];
rz(2.8475193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
