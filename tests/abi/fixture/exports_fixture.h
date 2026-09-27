/* Fixture for abi_exports_fires: never compiled. One function the library exports, and one it
 * does not. */
VCAL_API uint32_t VCAL_CALL vcal_abi_version(void);
VCAL_API vcal_status VCAL_CALL vcal_not_exported(void);
