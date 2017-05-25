
package us.kbase.kbstringtie;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: StringTieInput</p>
 * <pre>
 * required params:
 * assembly_ref: Alignment object reference
 * expression_set_name: ExpressionSet object name and output file header
 * workspace_name: the name of the workspace it gets saved to.
 * optional params:
 * num_threads: number of processing threads
 * junction_base: junctions that don't have spliced reads
 * junction_coverage: junction coverage
 * disable_trimming: disables trimming at the ends of the assembled transcripts
 * min_locus_gap_sep_value: minimum locus gap separation value
 * ballgown_mode: enables the output of Ballgown input table files
 * skip_reads_with_no_ref: reads with no reference will be skipped
 * maximum_fraction: maximum fraction of muliple-location-mapped reads
 * label: prefix for the name of the output transcripts
 * min_length: minimum length allowed for the predicted transcripts
 * min_read_coverage: minimum input transcript coverage
 * min_isoform_abundance: minimum isoform abundance
 * merge: set transcript merge mode
 * ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "alignment_ref",
    "expression_object_name",
    "workspace_name",
    "merge",
    "num_threads",
    "junction_base",
    "junction_coverage",
    "disable_trimming",
    "min_locus_gap_sep_value",
    "ballgown_mode",
    "skip_reads_with_no_ref",
    "maximum_fraction",
    "label",
    "min_length",
    "min_read_coverage",
    "min_isoform_abundance"
})
public class StringTieInput {

    @JsonProperty("alignment_ref")
    private String alignmentRef;
    @JsonProperty("expression_object_name")
    private String expressionObjectName;
    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("merge")
    private Long merge;
    @JsonProperty("num_threads")
    private Long numThreads;
    @JsonProperty("junction_base")
    private Long junctionBase;
    @JsonProperty("junction_coverage")
    private Double junctionCoverage;
    @JsonProperty("disable_trimming")
    private Long disableTrimming;
    @JsonProperty("min_locus_gap_sep_value")
    private Long minLocusGapSepValue;
    @JsonProperty("ballgown_mode")
    private Long ballgownMode;
    @JsonProperty("skip_reads_with_no_ref")
    private Long skipReadsWithNoRef;
    @JsonProperty("maximum_fraction")
    private Double maximumFraction;
    @JsonProperty("label")
    private String label;
    @JsonProperty("min_length")
    private Long minLength;
    @JsonProperty("min_read_coverage")
    private Double minReadCoverage;
    @JsonProperty("min_isoform_abundance")
    private Double minIsoformAbundance;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("alignment_ref")
    public String getAlignmentRef() {
        return alignmentRef;
    }

    @JsonProperty("alignment_ref")
    public void setAlignmentRef(String alignmentRef) {
        this.alignmentRef = alignmentRef;
    }

    public StringTieInput withAlignmentRef(String alignmentRef) {
        this.alignmentRef = alignmentRef;
        return this;
    }

    @JsonProperty("expression_object_name")
    public String getExpressionObjectName() {
        return expressionObjectName;
    }

    @JsonProperty("expression_object_name")
    public void setExpressionObjectName(String expressionObjectName) {
        this.expressionObjectName = expressionObjectName;
    }

    public StringTieInput withExpressionObjectName(String expressionObjectName) {
        this.expressionObjectName = expressionObjectName;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public StringTieInput withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("merge")
    public Long getMerge() {
        return merge;
    }

    @JsonProperty("merge")
    public void setMerge(Long merge) {
        this.merge = merge;
    }

    public StringTieInput withMerge(Long merge) {
        this.merge = merge;
        return this;
    }

    @JsonProperty("num_threads")
    public Long getNumThreads() {
        return numThreads;
    }

    @JsonProperty("num_threads")
    public void setNumThreads(Long numThreads) {
        this.numThreads = numThreads;
    }

    public StringTieInput withNumThreads(Long numThreads) {
        this.numThreads = numThreads;
        return this;
    }

    @JsonProperty("junction_base")
    public Long getJunctionBase() {
        return junctionBase;
    }

    @JsonProperty("junction_base")
    public void setJunctionBase(Long junctionBase) {
        this.junctionBase = junctionBase;
    }

    public StringTieInput withJunctionBase(Long junctionBase) {
        this.junctionBase = junctionBase;
        return this;
    }

    @JsonProperty("junction_coverage")
    public Double getJunctionCoverage() {
        return junctionCoverage;
    }

    @JsonProperty("junction_coverage")
    public void setJunctionCoverage(Double junctionCoverage) {
        this.junctionCoverage = junctionCoverage;
    }

    public StringTieInput withJunctionCoverage(Double junctionCoverage) {
        this.junctionCoverage = junctionCoverage;
        return this;
    }

    @JsonProperty("disable_trimming")
    public Long getDisableTrimming() {
        return disableTrimming;
    }

    @JsonProperty("disable_trimming")
    public void setDisableTrimming(Long disableTrimming) {
        this.disableTrimming = disableTrimming;
    }

    public StringTieInput withDisableTrimming(Long disableTrimming) {
        this.disableTrimming = disableTrimming;
        return this;
    }

    @JsonProperty("min_locus_gap_sep_value")
    public Long getMinLocusGapSepValue() {
        return minLocusGapSepValue;
    }

    @JsonProperty("min_locus_gap_sep_value")
    public void setMinLocusGapSepValue(Long minLocusGapSepValue) {
        this.minLocusGapSepValue = minLocusGapSepValue;
    }

    public StringTieInput withMinLocusGapSepValue(Long minLocusGapSepValue) {
        this.minLocusGapSepValue = minLocusGapSepValue;
        return this;
    }

    @JsonProperty("ballgown_mode")
    public Long getBallgownMode() {
        return ballgownMode;
    }

    @JsonProperty("ballgown_mode")
    public void setBallgownMode(Long ballgownMode) {
        this.ballgownMode = ballgownMode;
    }

    public StringTieInput withBallgownMode(Long ballgownMode) {
        this.ballgownMode = ballgownMode;
        return this;
    }

    @JsonProperty("skip_reads_with_no_ref")
    public Long getSkipReadsWithNoRef() {
        return skipReadsWithNoRef;
    }

    @JsonProperty("skip_reads_with_no_ref")
    public void setSkipReadsWithNoRef(Long skipReadsWithNoRef) {
        this.skipReadsWithNoRef = skipReadsWithNoRef;
    }

    public StringTieInput withSkipReadsWithNoRef(Long skipReadsWithNoRef) {
        this.skipReadsWithNoRef = skipReadsWithNoRef;
        return this;
    }

    @JsonProperty("maximum_fraction")
    public Double getMaximumFraction() {
        return maximumFraction;
    }

    @JsonProperty("maximum_fraction")
    public void setMaximumFraction(Double maximumFraction) {
        this.maximumFraction = maximumFraction;
    }

    public StringTieInput withMaximumFraction(Double maximumFraction) {
        this.maximumFraction = maximumFraction;
        return this;
    }

    @JsonProperty("label")
    public String getLabel() {
        return label;
    }

    @JsonProperty("label")
    public void setLabel(String label) {
        this.label = label;
    }

    public StringTieInput withLabel(String label) {
        this.label = label;
        return this;
    }

    @JsonProperty("min_length")
    public Long getMinLength() {
        return minLength;
    }

    @JsonProperty("min_length")
    public void setMinLength(Long minLength) {
        this.minLength = minLength;
    }

    public StringTieInput withMinLength(Long minLength) {
        this.minLength = minLength;
        return this;
    }

    @JsonProperty("min_read_coverage")
    public Double getMinReadCoverage() {
        return minReadCoverage;
    }

    @JsonProperty("min_read_coverage")
    public void setMinReadCoverage(Double minReadCoverage) {
        this.minReadCoverage = minReadCoverage;
    }

    public StringTieInput withMinReadCoverage(Double minReadCoverage) {
        this.minReadCoverage = minReadCoverage;
        return this;
    }

    @JsonProperty("min_isoform_abundance")
    public Double getMinIsoformAbundance() {
        return minIsoformAbundance;
    }

    @JsonProperty("min_isoform_abundance")
    public void setMinIsoformAbundance(Double minIsoformAbundance) {
        this.minIsoformAbundance = minIsoformAbundance;
    }

    public StringTieInput withMinIsoformAbundance(Double minIsoformAbundance) {
        this.minIsoformAbundance = minIsoformAbundance;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((((((((((((("StringTieInput"+" [alignmentRef=")+ alignmentRef)+", expressionObjectName=")+ expressionObjectName)+", workspaceName=")+ workspaceName)+", merge=")+ merge)+", numThreads=")+ numThreads)+", junctionBase=")+ junctionBase)+", junctionCoverage=")+ junctionCoverage)+", disableTrimming=")+ disableTrimming)+", minLocusGapSepValue=")+ minLocusGapSepValue)+", ballgownMode=")+ ballgownMode)+", skipReadsWithNoRef=")+ skipReadsWithNoRef)+", maximumFraction=")+ maximumFraction)+", label=")+ label)+", minLength=")+ minLength)+", minReadCoverage=")+ minReadCoverage)+", minIsoformAbundance=")+ minIsoformAbundance)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
