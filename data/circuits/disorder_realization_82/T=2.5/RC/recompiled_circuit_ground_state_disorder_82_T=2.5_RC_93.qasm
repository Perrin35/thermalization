OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3759484) q[0];
sx q[0];
rz(3.5557616) q[0];
sx q[0];
rz(10.226305) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(2.2771775) q[1];
sx q[1];
rz(9.519001) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163001) q[0];
sx q[0];
rz(-2.2803118) q[0];
sx q[0];
rz(0.50523357) q[0];
rz(-pi) q[1];
rz(1.6802189) q[2];
sx q[2];
rz(-2.8673807) q[2];
sx q[2];
rz(0.69137979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6871145) q[1];
sx q[1];
rz(-1.2839437) q[1];
sx q[1];
rz(-1.0832562) q[1];
rz(-pi) q[2];
rz(2.9889936) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(1.2945557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6039383) q[2];
sx q[2];
rz(-0.64017576) q[2];
sx q[2];
rz(1.630265) q[2];
rz(-0.18621914) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(2.490624) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52361012) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(0.13667662) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-2.3725841) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5154607) q[0];
sx q[0];
rz(-0.40053408) q[0];
sx q[0];
rz(-0.35450165) q[0];
rz(-2.7203975) q[2];
sx q[2];
rz(-2.4026818) q[2];
sx q[2];
rz(0.74904672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8254346) q[1];
sx q[1];
rz(-2.3806751) q[1];
sx q[1];
rz(-2.9050211) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8212998) q[3];
sx q[3];
rz(-0.13938306) q[3];
sx q[3];
rz(-1.936862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-2.365999) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8785777) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(2.0468792) q[0];
rz(1.8393983) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(-0.0040231752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4132679) q[0];
sx q[0];
rz(-2.2906429) q[0];
sx q[0];
rz(-2.9664449) q[0];
rz(-pi) q[1];
rz(0.74538576) q[2];
sx q[2];
rz(-2.2507387) q[2];
sx q[2];
rz(0.58733515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6371256) q[1];
sx q[1];
rz(-1.5404623) q[1];
sx q[1];
rz(-2.8843969) q[1];
rz(-pi) q[2];
rz(-1.8173006) q[3];
sx q[3];
rz(-0.87216264) q[3];
sx q[3];
rz(-2.7072631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5912938) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(2.6934521) q[2];
rz(0.48629931) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(-1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(1.9122596) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(-2.5965447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16825039) q[0];
sx q[0];
rz(-2.0297883) q[0];
sx q[0];
rz(-0.58990546) q[0];
x q[1];
rz(-2.2322234) q[2];
sx q[2];
rz(-1.0535568) q[2];
sx q[2];
rz(1.478372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74074827) q[1];
sx q[1];
rz(-1.6312338) q[1];
sx q[1];
rz(-1.7857741) q[1];
x q[2];
rz(-2.0470951) q[3];
sx q[3];
rz(-2.3102488) q[3];
sx q[3];
rz(0.0090713105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58686078) q[2];
sx q[2];
rz(-1.502159) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(1.933291) q[3];
sx q[3];
rz(-2.8390563) q[3];
sx q[3];
rz(1.6944073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(0.30715352) q[1];
sx q[1];
rz(-0.75633621) q[1];
sx q[1];
rz(-0.53877962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2870432) q[0];
sx q[0];
rz(-1.3417305) q[0];
sx q[0];
rz(1.7160709) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42947763) q[2];
sx q[2];
rz(-1.4444077) q[2];
sx q[2];
rz(-0.9867368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9591892) q[1];
sx q[1];
rz(-0.99970523) q[1];
sx q[1];
rz(1.1491999) q[1];
rz(-pi) q[2];
rz(2.6461824) q[3];
sx q[3];
rz(-0.79548478) q[3];
sx q[3];
rz(-2.0919679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6946081) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(1.2896607) q[2];
rz(-1.5223711) q[3];
sx q[3];
rz(-2.3963942) q[3];
sx q[3];
rz(1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(-2.4094021) q[1];
sx q[1];
rz(-0.79473549) q[1];
sx q[1];
rz(-2.5644459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2271778) q[0];
sx q[0];
rz(-3.1127226) q[0];
sx q[0];
rz(-1.6976852) q[0];
x q[1];
rz(0.24546282) q[2];
sx q[2];
rz(-2.5405014) q[2];
sx q[2];
rz(-2.9562151) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0476372) q[1];
sx q[1];
rz(-3.006662) q[1];
sx q[1];
rz(2.5173902) q[1];
rz(-0.34586819) q[3];
sx q[3];
rz(-2.3551999) q[3];
sx q[3];
rz(-2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.131549) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.6183759) q[2];
rz(-3.0736382) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(-2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0129358) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(2.3701684) q[0];
rz(-0.5386638) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(-2.0753986) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5504042) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(-1.2729086) q[0];
x q[1];
rz(1.9692405) q[2];
sx q[2];
rz(-0.5486998) q[2];
sx q[2];
rz(-1.0469701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2354338) q[1];
sx q[1];
rz(-1.3536436) q[1];
sx q[1];
rz(1.557334) q[1];
x q[2];
rz(-2.8071219) q[3];
sx q[3];
rz(-1.3190509) q[3];
sx q[3];
rz(-1.2125963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44925877) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(-0.44245693) q[2];
rz(0.19065204) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.8316487) q[0];
rz(0.39778057) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(1.3936874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3041954) q[0];
sx q[0];
rz(-2.1050354) q[0];
sx q[0];
rz(-2.6287352) q[0];
rz(2.2289235) q[2];
sx q[2];
rz(-2.7670112) q[2];
sx q[2];
rz(-0.41824579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19119769) q[1];
sx q[1];
rz(-2.0346332) q[1];
sx q[1];
rz(1.9056803) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3721136) q[3];
sx q[3];
rz(-1.4666838) q[3];
sx q[3];
rz(-2.7840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.028367793) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(0.69463378) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46591127) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(-2.1821816) q[0];
rz(0.77964669) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(0.27228212) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36764964) q[0];
sx q[0];
rz(-3.127357) q[0];
sx q[0];
rz(1.9444501) q[0];
rz(1.2113038) q[2];
sx q[2];
rz(-0.73795107) q[2];
sx q[2];
rz(-2.8518845) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6710558) q[1];
sx q[1];
rz(-1.1644286) q[1];
sx q[1];
rz(-1.0419091) q[1];
x q[2];
rz(-2.394478) q[3];
sx q[3];
rz(-2.0339587) q[3];
sx q[3];
rz(2.2833191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.73725629) q[2];
sx q[2];
rz(-1.6360838) q[2];
sx q[2];
rz(0.20757248) q[2];
rz(3.0370144) q[3];
sx q[3];
rz(-2.6360377) q[3];
sx q[3];
rz(1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.680147) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(1.0567868) q[0];
rz(-2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(-2.812885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8037107) q[0];
sx q[0];
rz(-2.0265371) q[0];
sx q[0];
rz(-2.7720939) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6255837) q[2];
sx q[2];
rz(-1.6431) q[2];
sx q[2];
rz(2.876578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.199177) q[1];
sx q[1];
rz(-2.7991382) q[1];
sx q[1];
rz(-1.5210995) q[1];
x q[2];
rz(2.8601951) q[3];
sx q[3];
rz(-0.79165047) q[3];
sx q[3];
rz(-2.7782914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1591961) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(2.0170434) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(-0.92397773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9818253) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(-0.81623296) q[1];
sx q[1];
rz(-1.7385794) q[1];
sx q[1];
rz(2.0508918) q[1];
rz(-2.6098567) q[2];
sx q[2];
rz(-2.3219863) q[2];
sx q[2];
rz(1.815997) q[2];
rz(0.40533752) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
