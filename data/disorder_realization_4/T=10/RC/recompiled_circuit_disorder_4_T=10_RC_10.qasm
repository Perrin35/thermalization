OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(-0.43570575) q[0];
sx q[0];
rz(-2.2154007) q[0];
rz(-1.1821094) q[1];
sx q[1];
rz(3.8745772) q[1];
sx q[1];
rz(12.193845) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46967888) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(-1.8125305) q[0];
x q[1];
rz(-0.1442285) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(0.47338212) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8718611) q[1];
sx q[1];
rz(-1.4565804) q[1];
sx q[1];
rz(1.7608587) q[1];
x q[2];
rz(3.0285809) q[3];
sx q[3];
rz(-1.7859308) q[3];
sx q[3];
rz(2.5290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7131876) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(-2.2170846) q[2];
rz(1.6690147) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(-1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(-1.5361319) q[0];
rz(-2.9470782) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(-3.0867192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957735) q[0];
sx q[0];
rz(-1.6841444) q[0];
sx q[0];
rz(0.093009526) q[0];
x q[1];
rz(1.0286249) q[2];
sx q[2];
rz(-1.0222058) q[2];
sx q[2];
rz(2.2749535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.730719) q[1];
sx q[1];
rz(-1.5667856) q[1];
sx q[1];
rz(-1.2282759) q[1];
rz(-pi) q[2];
rz(-1.7826471) q[3];
sx q[3];
rz(-2.2329674) q[3];
sx q[3];
rz(0.3609095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.6595718) q[2];
rz(-0.99059087) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822587) q[0];
sx q[0];
rz(-2.6222836) q[0];
sx q[0];
rz(-0.3749795) q[0];
rz(0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.8050271) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35555392) q[0];
sx q[0];
rz(-2.9125014) q[0];
sx q[0];
rz(-1.8285719) q[0];
rz(-2.5518718) q[2];
sx q[2];
rz(-1.9460742) q[2];
sx q[2];
rz(-0.42082618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3805483) q[1];
sx q[1];
rz(-2.3420463) q[1];
sx q[1];
rz(-1.9995081) q[1];
x q[2];
rz(-0.7671719) q[3];
sx q[3];
rz(-1.2734405) q[3];
sx q[3];
rz(-1.5945827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(-1.864795) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643726) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(1.9112446) q[1];
sx q[1];
rz(-2.3760445) q[1];
sx q[1];
rz(-1.4039325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0899635) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(-2.7882663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4594175) q[2];
sx q[2];
rz(-2.1795159) q[2];
sx q[2];
rz(-1.9145554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0701323) q[1];
sx q[1];
rz(-1.4590108) q[1];
sx q[1];
rz(-1.2791355) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0560533) q[3];
sx q[3];
rz(-2.2604058) q[3];
sx q[3];
rz(-0.52348189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4924865) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(-2.9283294) q[2];
rz(-3.1212741) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(-0.4030574) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(1.0158585) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(3.1075081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406697) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(-2.7435061) q[0];
x q[1];
rz(-1.3047406) q[2];
sx q[2];
rz(-2.1286466) q[2];
sx q[2];
rz(-1.7709874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2195697) q[1];
sx q[1];
rz(-1.2864188) q[1];
sx q[1];
rz(-3.1344096) q[1];
x q[2];
rz(-2.2250416) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.3548917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(-1.9704698) q[2];
rz(-1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(-3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(0.43584287) q[0];
rz(1.1054989) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(-1.2058535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642731) q[0];
sx q[0];
rz(-0.84007971) q[0];
sx q[0];
rz(0.3198091) q[0];
rz(1.7370991) q[2];
sx q[2];
rz(-0.096509343) q[2];
sx q[2];
rz(-1.6627179) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7534415) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(-0.053877342) q[1];
rz(-1.9704351) q[3];
sx q[3];
rz(-1.6013147) q[3];
sx q[3];
rz(0.37351028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-3.0899866) q[2];
rz(2.7339325) q[3];
sx q[3];
rz(-1.1838653) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62149858) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.4173123) q[1];
sx q[1];
rz(2.6046682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818031) q[0];
sx q[0];
rz(-1.4872695) q[0];
sx q[0];
rz(3.0755755) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9951018) q[2];
sx q[2];
rz(-1.8244201) q[2];
sx q[2];
rz(0.25073642) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9432224) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(-0.13065773) q[1];
x q[2];
rz(-1.7241931) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(3.1068387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9817104) q[2];
sx q[2];
rz(-1.8843702) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(0.28132176) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(2.8040335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71096703) q[0];
sx q[0];
rz(-1.0904878) q[0];
sx q[0];
rz(-1.7947012) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0701596) q[2];
sx q[2];
rz(-2.2087503) q[2];
sx q[2];
rz(0.86779867) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2355289) q[1];
sx q[1];
rz(-1.4396832) q[1];
sx q[1];
rz(0.99793418) q[1];
x q[2];
rz(1.6345331) q[3];
sx q[3];
rz(-0.81634854) q[3];
sx q[3];
rz(3.0626429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(-2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(2.3601941) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18701126) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(0.22098456) q[0];
rz(0.9221319) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(-2.1386713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10119469) q[0];
sx q[0];
rz(-1.4324491) q[0];
sx q[0];
rz(-0.25427108) q[0];
rz(2.6384764) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(-2.2803277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6332615) q[1];
sx q[1];
rz(-1.08053) q[1];
sx q[1];
rz(-2.0983178) q[1];
rz(-pi) q[2];
rz(-2.2569611) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(-2.6209944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(-2.6878099) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145988) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(-0.19113834) q[0];
rz(2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(0.89231649) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494209) q[0];
sx q[0];
rz(-2.3755709) q[0];
sx q[0];
rz(1.8269405) q[0];
rz(-pi) q[1];
rz(1.4453663) q[2];
sx q[2];
rz(-2.4969366) q[2];
sx q[2];
rz(0.21412011) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72496966) q[1];
sx q[1];
rz(-1.5048711) q[1];
sx q[1];
rz(-0.23857393) q[1];
rz(-pi) q[2];
rz(-2.7097706) q[3];
sx q[3];
rz(-1.8027657) q[3];
sx q[3];
rz(-2.1524129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(1.9178948) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(2.3971476) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214462) q[0];
sx q[0];
rz(-2.2059724) q[0];
sx q[0];
rz(-1.7511517) q[0];
rz(1.4795115) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(1.66172) q[2];
sx q[2];
rz(-0.29279136) q[2];
sx q[2];
rz(1.472483) q[2];
rz(-0.039924351) q[3];
sx q[3];
rz(-1.6394686) q[3];
sx q[3];
rz(-2.6185262) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
