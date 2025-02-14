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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(-1.1535147) q[0];
rz(1.310362) q[1];
sx q[1];
rz(5.789776) q[1];
sx q[1];
rz(8.9763666) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759244) q[0];
sx q[0];
rz(-2.2184555) q[0];
sx q[0];
rz(1.5041385) q[0];
x q[1];
rz(-1.469606) q[2];
sx q[2];
rz(-1.8253583) q[2];
sx q[2];
rz(2.2806666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4223215) q[1];
sx q[1];
rz(-1.0485253) q[1];
sx q[1];
rz(2.2838255) q[1];
x q[2];
rz(2.6403994) q[3];
sx q[3];
rz(-1.6128916) q[3];
sx q[3];
rz(2.834109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6260234) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-0.62421978) q[2];
rz(-2.7624687) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17495951) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(-2.1061184) q[0];
rz(1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(-0.48286352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.567822) q[0];
sx q[0];
rz(-1.616524) q[0];
sx q[0];
rz(-3.0450495) q[0];
rz(-0.16629433) q[2];
sx q[2];
rz(-1.4539084) q[2];
sx q[2];
rz(2.6720195) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16115621) q[1];
sx q[1];
rz(-2.709678) q[1];
sx q[1];
rz(1.5932139) q[1];
rz(-0.32093908) q[3];
sx q[3];
rz(-1.4934936) q[3];
sx q[3];
rz(0.86107218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1172993) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.8710322) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(0.89699927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-0.73905149) q[0];
rz(2.1801379) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-2.3846073) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56102806) q[0];
sx q[0];
rz(-1.8317199) q[0];
sx q[0];
rz(-0.83103128) q[0];
rz(-pi) q[1];
rz(-2.4503373) q[2];
sx q[2];
rz(-1.2757841) q[2];
sx q[2];
rz(-1.472109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0464335) q[1];
sx q[1];
rz(-0.94645247) q[1];
sx q[1];
rz(1.4051564) q[1];
rz(-pi) q[2];
rz(-1.2631031) q[3];
sx q[3];
rz(-1.6374) q[3];
sx q[3];
rz(0.25552017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-0.91705489) q[2];
sx q[2];
rz(-2.4364831) q[2];
rz(-0.32143587) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-2.7307935) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.131677) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(-1.2257082) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(3.0753678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787251) q[0];
sx q[0];
rz(-1.6604794) q[0];
sx q[0];
rz(1.8171726) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29860626) q[2];
sx q[2];
rz(-0.82393194) q[2];
sx q[2];
rz(0.9048942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5410616) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(-2.1670684) q[1];
x q[2];
rz(-0.36599593) q[3];
sx q[3];
rz(-1.4533236) q[3];
sx q[3];
rz(2.7155196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(1.0062086) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(-0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202268) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(1.5059858) q[0];
rz(1.6141363) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(-1.45586) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1560178) q[0];
sx q[0];
rz(-2.394372) q[0];
sx q[0];
rz(1.78563) q[0];
x q[1];
rz(-1.4541995) q[2];
sx q[2];
rz(-1.3449838) q[2];
sx q[2];
rz(1.9436398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29042915) q[1];
sx q[1];
rz(-2.4896087) q[1];
sx q[1];
rz(-0.033520582) q[1];
x q[2];
rz(2.5334356) q[3];
sx q[3];
rz(-1.4440315) q[3];
sx q[3];
rz(-0.34846445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(0.91147649) q[2];
rz(-2.2677926) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-3.0112596) q[0];
rz(0.48769543) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(0.74388751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2206534) q[0];
sx q[0];
rz(-1.6531117) q[0];
sx q[0];
rz(-3.0916075) q[0];
x q[1];
rz(-3.0024372) q[2];
sx q[2];
rz(-0.94991131) q[2];
sx q[2];
rz(-0.54229743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0186613) q[1];
sx q[1];
rz(-2.5187188) q[1];
sx q[1];
rz(-1.0708924) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3906995) q[3];
sx q[3];
rz(-2.0493747) q[3];
sx q[3];
rz(-0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1193715) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(-2.1807097) q[2];
rz(-2.7458701) q[3];
sx q[3];
rz(-1.3362249) q[3];
sx q[3];
rz(-3.0254288) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4294089) q[0];
sx q[0];
rz(-1.3954221) q[0];
sx q[0];
rz(-0.87923572) q[0];
rz(3.0839594) q[2];
sx q[2];
rz(-2.1890292) q[2];
sx q[2];
rz(2.1098441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8128374) q[1];
sx q[1];
rz(-2.8649674) q[1];
sx q[1];
rz(0.098621086) q[1];
rz(-pi) q[2];
rz(1.3052528) q[3];
sx q[3];
rz(-1.2528386) q[3];
sx q[3];
rz(1.4972655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0543694) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(1.8325904) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.8488309) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(-2.8344179) q[0];
rz(1.81709) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-0.90075341) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9884315) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(-0.26217006) q[0];
x q[1];
rz(-1.0020761) q[2];
sx q[2];
rz(-2.4042712) q[2];
sx q[2];
rz(-1.7676644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5430205) q[1];
sx q[1];
rz(-1.9024339) q[1];
sx q[1];
rz(0.33473067) q[1];
rz(-pi) q[2];
rz(-1.7124655) q[3];
sx q[3];
rz(-2.1974539) q[3];
sx q[3];
rz(3.0350445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8073392) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(2.0452512) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39356247) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(-0.50931859) q[0];
rz(-2.9604984) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(0.96955713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60107952) q[0];
sx q[0];
rz(-0.76427312) q[0];
sx q[0];
rz(-2.918052) q[0];
rz(-pi) q[1];
rz(-3.0907057) q[2];
sx q[2];
rz(-1.0641838) q[2];
sx q[2];
rz(1.8614872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6020376) q[1];
sx q[1];
rz(-1.4629211) q[1];
sx q[1];
rz(-2.8007568) q[1];
rz(-pi) q[2];
rz(-0.84253715) q[3];
sx q[3];
rz(-1.304783) q[3];
sx q[3];
rz(-2.0913948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-1.0815557) q[2];
rz(-2.9099416) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-0.20802465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.76081) q[0];
sx q[0];
rz(-0.2247227) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(-2.6864247) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(-0.761935) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2357016) q[0];
sx q[0];
rz(-2.3755223) q[0];
sx q[0];
rz(2.7188042) q[0];
rz(-0.81934118) q[2];
sx q[2];
rz(-1.0723253) q[2];
sx q[2];
rz(0.33453951) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40574351) q[1];
sx q[1];
rz(-2.4824576) q[1];
sx q[1];
rz(3.0856611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1033789) q[3];
sx q[3];
rz(-1.5776792) q[3];
sx q[3];
rz(-0.21871834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(-0.90395149) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91697964) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(1.0820816) q[2];
sx q[2];
rz(-2.6398224) q[2];
sx q[2];
rz(2.421438) q[2];
rz(-0.42999646) q[3];
sx q[3];
rz(-1.5734133) q[3];
sx q[3];
rz(-1.6071241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
